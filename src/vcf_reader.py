import logging
from typing import Dict, Any, Tuple, Optional

import allel

from base.reference_site import ReferenceSite
from base.util import strip_prefix
from config.tool_config import ToolConfig
from base.constants import REF_CALL_ANNOTATION_STRING
from base.filter import SimpleCallFilter
from base.gene_coordinate import GeneCoordinate
from base.reference_assembly import ReferenceAssembly
from call_data import SimpleCallData, SimpleCall
from config.panel import Panel


class VcfReader(object):
    CHROMOSOME_FIELD_NAME = "variants/CHROM"
    POSITION_FIELD_NAME = "variants/POS"
    RS_IDS_FIELD_NAME = "variants/ID"
    REF_ALLELE_FIELD_NAME = "variants/REF"
    ALT_ALLELE_FIELD_NAME = "variants/ALT"
    FILTER_FIELD_NAME = "variants/FILTER"
    ANNOTATION_FIELD_NAME = "variants/ANN"
    GENOTYPE_FIELD_NAME = "calldata/GT"
    SAMPLE_FIELD_NAME = "samples"

    RS_ID_SEPARATOR = ";"
    RS_ID_EMPTY_INDICATOR = "."
    CODING_VARIANT_ANNOTATION_PREFIX = "c."
    NON_CODING_VARIANT_ANNOTATION_PREFIX = "n."

    def get_call_data(self, tool_config: ToolConfig, panel: Panel) -> SimpleCallData:
        variants = self.__get_variants_from_vcf(tool_config.vcf_path)
        if variants is not None:
            return self.__get_call_data_from_variants(
                variants, panel, tool_config.sample_r_id, tool_config.vcf_reference_assembly
            )
        else:
            logging.warning("No variants found in vcf")
            return SimpleCallData(frozenset(), tool_config.vcf_reference_assembly)

    def __get_variants_from_vcf(self, vcf: str) -> Optional[Dict[str, Any]]:
        # variants is None precisely when the vcf file has no variants
        try:
            variants = allel.read_vcf(vcf, fields="*")
        except IOError:
            raise FileNotFoundError("File " + vcf + " not found or cannot be opened.")
        return variants

    def __get_call_data_from_variants(
        self, variants: Dict[str, Any], panel: Panel, sample_r_id: str, vcf_reference_assembly: ReferenceAssembly
    ) -> SimpleCallData:
        filtered_calls = set()

        total_variant_count = self.__get_variant_count(variants)
        logging.info(f"VCF calls: {total_variant_count}")
        for call_index in range(total_variant_count):
            if not self.__filter_is_pass(call_index, variants):
                # Ignore all calls with filter != PASS
                continue

            call = self.__get_call_from_variants(call_index, sample_r_id, variants)
            if panel.is_relevant_to_panel(call, vcf_reference_assembly):
                filtered_calls.add(call)

        logging.info(f"VCF calls QC-PASS and matching panel: {len(filtered_calls)}")

        return SimpleCallData(frozenset(filtered_calls), vcf_reference_assembly)

    def __filter_is_pass(self, call_index: int, variants: Dict[str, Any]) -> bool:
        return bool(variants[f"{self.FILTER_FIELD_NAME}_PASS"][call_index])

    def __get_call_from_variants(self, call_index: int, sample_r_id: str, variants: Dict[str, Any]) -> SimpleCall:
        chromosome = self.__get_chromosome_from_variants(call_index, variants)
        position = self.__get_position_from_variants(call_index, variants)
        gene_coordinate = GeneCoordinate(chromosome, position)

        reference_allele = self.__get_reference_allele_from_variants(call_index, variants)
        alleles = self.__get_called_alleles_from_variants(call_index, sample_r_id, variants)
        rs_ids = self.__get_rs_ids_from_variants(call_index, variants)

        if alleles == (reference_allele, reference_allele):
            gene_name, _ = self.__get_gene_name_and_variant_annotation_from_variants(call_index, variants)
            variant_annotation = REF_CALL_ANNOTATION_STRING
        else:
            gene_name, variant_annotation = self.__get_gene_name_and_variant_annotation_from_variants(
                call_index, variants
            )

        call = SimpleCall(
            ReferenceSite(gene_coordinate, reference_allele),
            alleles,
            gene_name,
            rs_ids,
            variant_annotation,
            SimpleCallFilter.PASS,
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

    def __get_gene_name_and_variant_annotation_from_variants(
            self, call_index: int, variants: Dict[str, Any]
    ) -> Tuple[str, str]:
        complete_annotation = str(variants[self.ANNOTATION_FIELD_NAME][call_index])
        if complete_annotation:
            gene_name = complete_annotation.split("|")[3]
            full_variant_annotation = complete_annotation.split("|")[9]
            if full_variant_annotation.startswith(self.CODING_VARIANT_ANNOTATION_PREFIX):
                variant_annotation = strip_prefix(full_variant_annotation, self.CODING_VARIANT_ANNOTATION_PREFIX)
            elif full_variant_annotation.startswith(self.NON_CODING_VARIANT_ANNOTATION_PREFIX):
                variant_annotation = strip_prefix(full_variant_annotation, self.NON_CODING_VARIANT_ANNOTATION_PREFIX)
            else:
                raise ValueError(f"Unexpected annotation prefix: {full_variant_annotation}")
        else:
            # TODO: solve this in a better way, or assert later that none of the returned simple calls have this
            gene_name = ""
            variant_annotation = ""

        return gene_name, variant_annotation

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
