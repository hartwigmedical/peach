import logging
from enum import Enum, auto
from typing import Dict, Any, Tuple, Optional

import allel

from util.reference_site import ReferenceSite
from tool_config import ToolConfig
from util.constants import REF_CALL_ANNOTATION_STRING
from util.filter import VcfCallFilter
from util.gene_coordinate import GeneCoordinate
from util.reference_assembly import ReferenceAssembly
from calls.vcf_call import VcfCall, VcfCallData
from panel.panel import Panel


class AnnotationType(Enum):
    # For calls with respect to the vcf reference assembly
    SNPEFF = auto()
    PAVE = auto()
    NONE = auto()


class VcfReader(object):
    CHROMOSOME_FIELD_NAME = "variants/CHROM"
    POSITION_FIELD_NAME = "variants/POS"
    RS_IDS_FIELD_NAME = "variants/ID"
    REF_ALLELE_FIELD_NAME = "variants/REF"
    ALT_ALLELE_FIELD_NAME = "variants/ALT"
    FILTER_FIELD_NAME = "variants/FILTER"
    GENOTYPE_FIELD_NAME = "calldata/GT"
    SAMPLE_FIELD_NAME = "samples"

    SNPEFF_ANNOTATION_FIELD_NAME = "variants/ANN"
    SNPEFF_ANNOTATION_SEPARATOR = "|"
    PAVE_ANNOTATION_FIELD_NAME = "variants/PAVE_TI"
    PAVE_ANNOTATION_SEPARATOR = "|"

    RS_ID_SEPARATOR = ";"
    RS_ID_EMPTY_INDICATOR = "."
    CODING_VARIANT_ANNOTATION_PREFIX = "c."
    NON_CODING_VARIANT_ANNOTATION_PREFIX = "n."

    def get_call_data(self, tool_config: ToolConfig, panel: Panel) -> VcfCallData:
        variants = self.__get_variants_from_vcf(tool_config.vcf_path)
        if variants is not None:
            return self.__get_call_data_from_variants(
                variants, panel, tool_config.sample_r_id, tool_config.vcf_reference_assembly
            )
        else:
            logging.warning("No variants found in vcf")
            return VcfCallData(frozenset(), tool_config.vcf_reference_assembly)

    def __get_variants_from_vcf(self, vcf: str) -> Optional[Dict[str, Any]]:
        # variants is None precisely when the VCF file has no variants
        try:
            variants = allel.read_vcf(vcf, fields="*", numbers={self.PAVE_ANNOTATION_FIELD_NAME: 1000})
        except IOError:
            raise FileNotFoundError(f"File {vcf} not found or cannot be opened.")
        return variants

    def __get_call_data_from_variants(
        self, variants: Dict[str, Any], panel: Panel, sample_r_id: str, vcf_reference_assembly: ReferenceAssembly
    ) -> VcfCallData:
        total_variant_count = self.__get_variant_count(variants)
        logging.info(f"VCF calls: {total_variant_count}")

        annotation_type = self.__get_annotation_type(variants)

        filtered_calls = set()
        for call_index in range(total_variant_count):
            if not self.__filter_is_pass(call_index, variants):
                # Ignore all calls with filter != PASS
                continue

            call = self.__get_call_from_variants(call_index, sample_r_id, variants, annotation_type, panel)
            if panel.is_relevant_to_panel(call, vcf_reference_assembly):
                filtered_calls.add(call)

        logging.info(f"VCF calls QC-PASS and matching panel: {len(filtered_calls)}")

        return VcfCallData(frozenset(filtered_calls), vcf_reference_assembly)

    def __filter_is_pass(self, call_index: int, variants: Dict[str, Any]) -> bool:
        return bool(variants[f"{self.FILTER_FIELD_NAME}_PASS"][call_index])

    def __get_call_from_variants(
            self,
            call_index: int,
            sample_r_id: str,
            variants: Dict[str, Any],
            annotation_type: AnnotationType,
            panel: Panel,
    ) -> VcfCall:
        chromosome = self.__get_chromosome_from_variants(call_index, variants)
        position = self.__get_position_from_variants(call_index, variants)
        gene_coordinate = GeneCoordinate(chromosome, position)

        reference_allele = self.__get_reference_allele_from_variants(call_index, variants)
        alleles = self.__get_called_alleles_from_variants(call_index, sample_r_id, variants)
        rs_ids = self.__get_rs_ids_from_variants(call_index, variants)

        variant_annotation: Optional[str]
        if alleles == (reference_allele, reference_allele):
            variant_annotation = REF_CALL_ANNOTATION_STRING
        else:
            variant_annotation = self.__get_variant_annotation_from_variants(
                call_index, variants, annotation_type, panel
            )

        call = VcfCall(
            ReferenceSite(gene_coordinate, reference_allele),
            alleles,
            rs_ids,
            variant_annotation,
            VcfCallFilter.PASS,
        )
        return call

    def __get_called_alleles_from_variants(
        self, call_index: int, sample_r_id: str, variants: Dict[str, Any]
    ) -> Tuple[str, str]:
        reference_allele = self.__get_reference_allele_from_variants(call_index, variants)
        sample_list = [str(sample_name) for sample_name in variants[self.SAMPLE_FIELD_NAME].tolist()]
        sample_index = sample_list.index(sample_r_id)
        genotype = [int(call) for call in variants[self.GENOTYPE_FIELD_NAME][call_index][sample_index].tolist()]
        alts = [str(allele) for allele in variants[self.ALT_ALLELE_FIELD_NAME][call_index]]
        if genotype == [0, 1]:
            alleles = (reference_allele, alts[0])
        elif genotype == [1, 1]:
            alleles = (alts[0], alts[0])
        elif genotype == [1, 2]:
            alleles = (alts[0], alts[1])
        elif genotype == [0, 0]:
            alleles = (reference_allele, reference_allele)
        else:
            error_msg = f"Genotype not found: {genotype}"
            raise ValueError(error_msg)
        return alleles

    def __get_variant_annotation_from_variants(
            self, call_index: int, variants: Dict[str, Any], annotation_type: AnnotationType, panel: Panel
    ) -> Optional[str]:
        full_variant_annotation: Optional[str]
        if annotation_type == AnnotationType.SNPEFF:
            complete_annotation = self.__get_snpeff_annotation_string(call_index, variants)
            if complete_annotation is not None:
                full_variant_annotation = complete_annotation.split(self.SNPEFF_ANNOTATION_SEPARATOR)[9]
            else:
                full_variant_annotation = None
        elif annotation_type == AnnotationType.PAVE:
            complete_annotation = self.__get_relevant_pave_annotation(call_index, variants, panel)
            if complete_annotation is not None:
                full_variant_annotation = complete_annotation.split(self.PAVE_ANNOTATION_SEPARATOR)[5]
            else:
                full_variant_annotation = None
        elif annotation_type == AnnotationType.NONE:
            complete_annotation = None
            full_variant_annotation = None
        else:
            raise ValueError(f"Unrecognized annotation type: {annotation_type}")

        variant_annotation: Optional[str]
        if full_variant_annotation is None or full_variant_annotation == "":
            variant_annotation = None
        elif full_variant_annotation.startswith(self.CODING_VARIANT_ANNOTATION_PREFIX):
            variant_annotation = self.__strip_prefix(full_variant_annotation, self.CODING_VARIANT_ANNOTATION_PREFIX)
        elif full_variant_annotation.startswith(self.NON_CODING_VARIANT_ANNOTATION_PREFIX):
            variant_annotation = self.__strip_prefix(full_variant_annotation, self.NON_CODING_VARIANT_ANNOTATION_PREFIX)
        else:
            error_msg = (
                f"Unexpected annotation prefix: "
                f"full_variant_annotation={full_variant_annotation}, "
                f"complete_annotation={complete_annotation}"
            )
            raise ValueError(error_msg)

        return variant_annotation

    def __get_relevant_pave_annotation(self, call_index: int, variants: Dict[str, Any], panel: Panel) -> Optional[str]:
        transcript_ids = panel.get_transcript_ids()
        if len(transcript_ids) != len(panel.get_genes()):
            error_msg = (
                f"For PEACH to handle PAVE annotations, the transcript ID needs to be set for all genes in the panel!"
            )
            raise ValueError(error_msg)
        transcript_id_to_annotation = self.__get_transcript_id_to_pave_annotation(call_index, variants)
        common_transcript_ids = transcript_ids.intersection(set(transcript_id_to_annotation.keys()))

        pave_annotation: Optional[str]
        if len(common_transcript_ids) == 1:
            pave_annotation = transcript_id_to_annotation[common_transcript_ids.pop()]
        elif not common_transcript_ids:
            pave_annotation = None
        else:
            error_msg = (
                f"Call has annotation for multiple transcripts of interest. "
                f"PEACH cannot handle genes with overlapping transcripts. "
                f"common_transcript_ids={common_transcript_ids}"
            )
            raise ValueError(error_msg)
        return pave_annotation

    def __get_transcript_id_to_pave_annotation(self, call_index: int, variants: Dict[str, Any]) -> Dict[str, str]:
        all_annotations = [str(annotation) for annotation in variants[self.PAVE_ANNOTATION_FIELD_NAME][call_index]]
        transcript_id_to_annotation: Dict[str, str] = {}
        for annotation in all_annotations:
            transcript_id = self.__get_transcript_id_from_pave_annotation(annotation)
            if transcript_id is not None:
                if transcript_id in transcript_id_to_annotation.keys():
                    error_msg = (
                        f"Call is annotated with a transcript ID multiple times: "
                        f"transcript_id={transcript_id}, call_index={call_index}"
                    )
                    raise ValueError(error_msg)
                transcript_id_to_annotation[transcript_id] = annotation
        return transcript_id_to_annotation

    def __get_transcript_id_from_pave_annotation(self, annotation: str) -> Optional[str]:
        # The list of PAVE annotations is padded with empty strings so all calls have the same number of annotations
        if annotation:
            return annotation.split(self.PAVE_ANNOTATION_SEPARATOR)[2]
        else:
            return None

    def __get_snpeff_annotation_string(self, call_index: int, variants: Dict[str, Any]) -> Optional[str]:
        complete_annotation = str(variants[self.SNPEFF_ANNOTATION_FIELD_NAME][call_index])
        if complete_annotation:
            return complete_annotation
        else:
            return None

    def __get_rs_ids_from_variants(self, call_index: int, variants: Dict[str, Any]) -> Tuple[str, ...]:
        rs_ids_string = str(variants[self.RS_IDS_FIELD_NAME][call_index])
        if self.RS_ID_SEPARATOR in rs_ids_string:
            return tuple(str(rs) for rs in rs_ids_string.split(self.RS_ID_SEPARATOR) if rs.startswith("rs"))
        elif rs_ids_string == self.RS_ID_EMPTY_INDICATOR:
            return tuple()
        else:
            return (str(rs_ids_string),)

    def __get_reference_allele_from_variants(self, call_index: int, variants: Dict[str, Any]) -> str:
        return str(variants[self.REF_ALLELE_FIELD_NAME][call_index])

    def __get_position_from_variants(self, call_index: int, variants: Dict[str, Any]) -> int:
        return int(variants[self.POSITION_FIELD_NAME][call_index])

    def __get_chromosome_from_variants(self, call_index: int, variants: Dict[str, Any]) -> str:
        return str(variants[self.CHROMOSOME_FIELD_NAME][call_index])

    def __get_variant_count(self, variants: Dict[str, Any]) -> int:
        return len(variants[self.RS_IDS_FIELD_NAME])

    def __get_annotation_type(self, variants: Dict[str, Any]) -> AnnotationType:
        has_snpeff_annotation = self.SNPEFF_ANNOTATION_FIELD_NAME in variants.keys()
        has_pave_annotation = self.PAVE_ANNOTATION_FIELD_NAME in variants.keys()
        if has_snpeff_annotation and has_pave_annotation:
            logging.warning(
                f"Both SNPEFF and PAVE annotation detected. When both are present, PAVE annotation is preferred."
            )
            return AnnotationType.PAVE
        elif has_snpeff_annotation and not has_pave_annotation:
            return AnnotationType.SNPEFF
        elif not has_snpeff_annotation and has_pave_annotation:
            return AnnotationType.PAVE
        else:
            error_msg = f"No annotation detected in the VCF. Annotate the VCF with either SNPEFF or PAVE."
            raise ValueError(error_msg)
            # logging.warning(f"No annotation detected in the VCF. Annotation with SNPEFF or PAVE is preferred.")
            # return AnnotationType.NONE

    def __strip_prefix(self, string: str, prefix: str) -> str:
        if string.startswith(prefix):
            return string[len(prefix):]
        else:
            error_msg = (
                f"String does not start with expected prefix: string={string}, prefix={prefix}"
            )
            raise ValueError(error_msg)
