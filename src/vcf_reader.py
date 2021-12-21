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

    @classmethod
    def get_call_data(cls, tool_config: ToolConfig, panel: Panel) -> SimpleCallData:
        variants = cls.__get_variants_from_vcf(tool_config.vcf_path)
        if variants is not None:
            return cls.__get_call_data_from_variants(
                variants, panel, tool_config.sample_r_id, tool_config.vcf_reference_assembly
            )
        else:
            logging.warning("No variants found in vcf")
            return SimpleCallData(frozenset(), tool_config.vcf_reference_assembly)

    @classmethod
    def __get_variants_from_vcf(cls, vcf: str) -> Optional[Dict[str, Any]]:
        # variants is None precisely when the vcf file has no variants
        try:
            variants = allel.read_vcf(vcf, fields="*")
        except IOError:
            raise FileNotFoundError("File " + vcf + " not found or cannot be opened.")
        return variants

    @classmethod
    def __get_call_data_from_variants(
        cls, variants: Dict[str, Any], panel: Panel, sample_r_id: str, vcf_reference_assembly: ReferenceAssembly
    ) -> SimpleCallData:
        filtered_calls = set()

        total_variant_count = cls.__get_variant_count(variants)
        logging.info(f"VCF calls: {total_variant_count}")
        for call_index in range(total_variant_count):
            if not cls.__filter_is_pass(call_index, variants):
                # Ignore all calls with filter != PASS
                continue

            call = cls.__get_call_from_variants(call_index, sample_r_id, variants)
            if panel.get_relevant_rs_ids(call, vcf_reference_assembly):
                filtered_calls.add(call)

        logging.info(f"VCF calls QC-PASS and matching panel: {len(filtered_calls)}")

        return SimpleCallData(frozenset(filtered_calls), vcf_reference_assembly)

    @classmethod
    def __filter_is_pass(cls, call_index: int, variants: Dict[str, Any]) -> bool:
        return bool(variants[f"{cls.FILTER_FIELD_NAME}_PASS"][call_index])

    @classmethod
    def __get_call_from_variants(cls, call_index: int, sample_r_id: str, variants: Dict[str, Any]) -> SimpleCall:
        chromosome = cls.__get_chromosome_from_variants(call_index, variants)
        position = cls.__get_position_from_variants(call_index, variants)
        gene_coordinate = GeneCoordinate(chromosome, position)

        reference_allele = cls.__get_reference_allele_from_variants(call_index, variants)
        alleles = cls.__get_called_alleles_from_variants(call_index, sample_r_id, variants)
        rs_ids = cls.__get_rs_ids_from_variants(call_index, variants)

        if alleles == (reference_allele, reference_allele):
            gene_name, _ = cls.__get_gene_name_and_variant_annotation_from_variants(call_index, variants)
            variant_annotation = REF_CALL_ANNOTATION_STRING
        else:
            gene_name, variant_annotation = cls.__get_gene_name_and_variant_annotation_from_variants(
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

    @classmethod
    def __get_called_alleles_from_variants(
        cls, call_index: int, sample_r_id: str, variants: Dict[str, Any]
    ) -> Tuple[str, str]:
        reference_allele = cls.__get_reference_allele_from_variants(call_index, variants)
        sample_list = [str(sample_name) for sample_name in variants[cls.SAMPLE_FIELD_NAME].tolist()]
        sample_index = sample_list.index(sample_r_id)
        genotype = [int(call) for call in variants[cls.GENOTYPE_FIELD_NAME][call_index][sample_index].tolist()]
        alts = [str(allele) for allele in variants[cls.ALT_ALLELE_FIELD_NAME][call_index]]
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

    @classmethod
    def __get_gene_name_and_variant_annotation_from_variants(
            cls, call_index: int, variants: Dict[str, Any]
    ) -> Tuple[str, str]:
        complete_annotation = str(variants[cls.ANNOTATION_FIELD_NAME][call_index])
        if complete_annotation:
            gene_name = complete_annotation.split("|")[3]
            full_variant_annotation = complete_annotation.split("|")[9]
            if full_variant_annotation.startswith(cls.CODING_VARIANT_ANNOTATION_PREFIX):
                variant_annotation = strip_prefix(full_variant_annotation, cls.CODING_VARIANT_ANNOTATION_PREFIX)
            elif full_variant_annotation.startswith(cls.NON_CODING_VARIANT_ANNOTATION_PREFIX):
                variant_annotation = strip_prefix(full_variant_annotation, cls.NON_CODING_VARIANT_ANNOTATION_PREFIX)
            else:
                raise ValueError(f"Unexpected annotation prefix: {full_variant_annotation}")
        else:
            # TODO: solve this in a better way, or assert later that none of the returned simple calls have this
            gene_name = ""
            variant_annotation = ""

        return gene_name, variant_annotation

    @classmethod
    def __get_rs_ids_from_variants(cls, call_index: int, variants: Dict[str, Any]) -> Tuple[str, ...]:
        rs_ids_string = str(variants[cls.RS_IDS_FIELD_NAME][call_index])
        if cls.RS_ID_SEPARATOR in rs_ids_string:
            return tuple(str(rs) for rs in rs_ids_string.split(cls.RS_ID_SEPARATOR) if rs.startswith("rs"))
        elif rs_ids_string == cls.RS_ID_EMPTY_INDICATOR:
            return tuple()
        else:
            return (str(rs_ids_string),)

    @classmethod
    def __get_reference_allele_from_variants(cls, call_index: int, variants: Dict[str, Any]) -> str:
        return str(variants[cls.REF_ALLELE_FIELD_NAME][call_index])

    @classmethod
    def __get_position_from_variants(cls, call_index: int, variants: Dict[str, Any]) -> int:
        return int(variants[cls.POSITION_FIELD_NAME][call_index])

    @classmethod
    def __get_chromosome_from_variants(cls, call_index: int, variants: Dict[str, Any]) -> str:
        return str(variants[cls.CHROMOSOME_FIELD_NAME][call_index])

    @classmethod
    def __get_variant_count(cls, variants: Dict[str, Any]) -> int:
        return len(variants[cls.RS_IDS_FIELD_NAME])
