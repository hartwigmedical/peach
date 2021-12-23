import itertools
import logging
from typing import List, Dict, Collection, FrozenSet, Optional, Set

from base.constants import NORMAL_FUNCTION_STRING
from base.reference_assembly import ReferenceAssembly
from base.util import get_key_to_multiple_values
from panel.annotation import Annotation
from panel.drug_info import DrugInfo, assert_no_overlap_drug_names
from panel.haplotype import (
    Haplotype,
    assert_no_overlap_haplotype_names,
    assert_no_overlap_haplotype_variant_combinations,
)
from panel.rs_id_info import RsIdInfo


class GeneInfo(object):
    """This object is meant to be immutable"""

    def __init__(
        self,
        gene: str,
        wild_type_haplotype_name: str,
        transcript_id: Optional[str],
        haplotypes: FrozenSet[Haplotype],
        rs_id_infos: FrozenSet[RsIdInfo],
        drugs: FrozenSet[DrugInfo],
        rs_id_to_ref_seq_difference_annotation: Dict[str, Annotation],
    ) -> None:
        assert_no_overlap_haplotype_names(haplotypes, f"gene info for {gene}")
        assert_no_overlap_haplotype_variant_combinations(haplotypes, f"gene info for {gene}")
        assert_no_overlap_drug_names(drugs, f"GeneInfo json for {gene}")
        self.__assert_rs_ids_all_different(rs_id_infos)
        self.__assert_rs_id_infos_compatible(rs_id_infos)
        self.__assert_rs_id_infos_match_chromosome(rs_id_infos, gene)
        self.__assert_info_exists_for_all_rs_ids_in_haplotypes(haplotypes, rs_id_infos)
        self.__assert_variants_in_haplotypes_compatible_with_rs_id_infos(haplotypes, rs_id_infos)
        self.__assert_rs_ids_with_ref_seq_differences_match_annotations(
            rs_id_infos, rs_id_to_ref_seq_difference_annotation
        )

        if not haplotypes:
            logging.warning(f"No alternate haplotypes configured for gene {gene}\n")

        self.__gene = gene
        self.__wild_type_haplotype_name = wild_type_haplotype_name
        self.__transcript_id = transcript_id
        self.__haplotypes = haplotypes
        self.__rs_id_infos = rs_id_infos
        self.__drugs = drugs
        self.__rs_id_to_ref_seq_difference_annotation = rs_id_to_ref_seq_difference_annotation

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, GeneInfo)
            and self.__gene == other.__gene
            and self.__wild_type_haplotype_name == other.__wild_type_haplotype_name
            and self.__transcript_id == other.__transcript_id
            and self.__haplotypes == other.__haplotypes
            and self.__rs_id_infos == other.__rs_id_infos
            and self.__drugs == other.__drugs
            and self.__rs_id_to_ref_seq_difference_annotation == other.__rs_id_to_ref_seq_difference_annotation
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.__gene,
                self.__wild_type_haplotype_name,
                self.__transcript_id,
                self.__haplotypes,
                self.__rs_id_infos,
                self.__drugs,
                tuple(sorted(self.__rs_id_to_ref_seq_difference_annotation.items())),
            )
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"GeneInfo("
            f"gene={self.__gene!r}, "
            f"wild_type_haplotype_name={self.__wild_type_haplotype_name!r}, "
            f"transcript_id={self.__transcript_id!r}, "
            f"haplotypes={self.__haplotypes!r}, "
            f"rs_id_infos={self.__rs_id_infos!r}, "
            f"drugs={self.__drugs!r}, "
            f"rs_id_to_ref_seq_difference_annotation={self.__rs_id_to_ref_seq_difference_annotation!r}, "
            f")"
        )

    @property
    def gene(self) -> str:
        return self.__gene

    @property
    def wild_type_haplotype_name(self) -> str:
        return self.__wild_type_haplotype_name

    @property
    def transcript_id(self) -> Optional[str]:
        return self.__transcript_id

    @property
    def haplotypes(self) -> FrozenSet[Haplotype]:
        return self.__haplotypes

    @property
    def rs_id_infos(self) -> FrozenSet[RsIdInfo]:
        return self.__rs_id_infos

    @property
    def drugs(self) -> FrozenSet[DrugInfo]:
        return self.__drugs

    def has_ref_sequence_difference_annotation(self, rs_id: str) -> bool:
        return rs_id in self.__rs_id_to_ref_seq_difference_annotation.keys()

    def get_ref_sequence_difference_annotation(self, rs_id: str, reference_assembly: ReferenceAssembly) -> str:
        if reference_assembly == ReferenceAssembly.V37:
            return self.__rs_id_to_ref_seq_difference_annotation[rs_id].annotation_v37
        elif reference_assembly == ReferenceAssembly.V38:
            return self.__rs_id_to_ref_seq_difference_annotation[rs_id].annotation_v38
        else:
            raise NotImplementedError(f"Unrecognized reference assembly: {reference_assembly}")

    def get_haplotype_function(self, haplotype_name: str) -> str:
        if haplotype_name == self.__wild_type_haplotype_name:
            return NORMAL_FUNCTION_STRING
        else:
            return self.__get_haplotype(haplotype_name).function

    def get_prescription_url(self, drug_name: str) -> str:
        matching_drug_infos = [drug_info for drug_info in self.__drugs if drug_info.name == drug_name]
        if len(matching_drug_infos) == 1:
            return matching_drug_infos[0].url_prescription_info
        else:
            raise ValueError(f"Not exactly one matching drug info in panel for gene={self.gene}: drug_name={drug_name}")

    def get_rs_ids(self) -> Set[str]:
        return {rs_id_info.rs_id for rs_id_info in self.rs_id_infos}

    def __get_haplotype(self, haplotype_name: str) -> Haplotype:
        matching_haplotypes = [haplotype for haplotype in self.__haplotypes if haplotype.name == haplotype_name]
        assert len(matching_haplotypes) == 1, (
            f"No unique haplotype with name: {haplotype_name} for gene {self.gene}:\n"
            f"matching_haplotypes={matching_haplotypes}"
        )
        return matching_haplotypes[0]

    @staticmethod
    def __assert_info_exists_for_all_rs_ids_in_haplotypes(
        haplotypes: FrozenSet[Haplotype], rs_id_infos: FrozenSet[RsIdInfo]
    ) -> None:
        rs_ids_in_haplotypes = {variant.rs_id for haplotype in haplotypes for variant in haplotype.variants}
        rs_ids_with_info = {info.rs_id for info in rs_id_infos}
        if not rs_ids_in_haplotypes.issubset(rs_ids_with_info):
            rs_ids_without_info = rs_ids_in_haplotypes.difference(rs_ids_with_info)
            error_msg = f"No info available for some of the rs ids in the haplotypes. Rs ids: {rs_ids_without_info}"
            raise ValueError(error_msg)

    @staticmethod
    def __assert_rs_ids_with_ref_seq_differences_match_annotations(
        rs_id_infos: FrozenSet[RsIdInfo], rs_id_to_ref_seq_difference_annotation: Dict[str, Annotation]
    ) -> None:
        rs_ids_from_infos = {info.rs_id for info in rs_id_infos if info.has_reference_sequence_difference()}
        rs_ids_from_annotation = set(rs_id_to_ref_seq_difference_annotation.keys())
        if rs_ids_from_infos != rs_ids_from_annotation:
            rs_ids_with_only_info = rs_ids_from_infos.difference(rs_ids_from_annotation)
            rs_ids_with_only_annotation = rs_ids_from_annotation.difference(rs_ids_from_infos)
            error_msg = (
                f"Rs ids with differences between v37 and v38 do not match "
                f"the rs ids with annotations for these differences. "
                f"Only info: {rs_ids_with_only_info} "
                f"Only annotation: {rs_ids_with_only_annotation}"
            )
            raise ValueError(error_msg)

    @staticmethod
    def __assert_rs_id_infos_compatible(rs_id_infos: FrozenSet[RsIdInfo]) -> None:
        for left, right in itertools.combinations(rs_id_infos, 2):
            if not left.is_compatible(right):
                error_msg = f"Incompatible rs id infos in gene info. left: {left}, right: {right}"
                raise ValueError(error_msg)

    @staticmethod
    def __assert_rs_id_infos_match_chromosome(rs_id_infos: FrozenSet[RsIdInfo], gene: str) -> None:
        v37_chromosomes = {info.reference_site_v37.start_coordinate.chromosome for info in rs_id_infos}
        v38_chromosomes = {info.reference_site_v38.start_coordinate.chromosome for info in rs_id_infos}
        if len(v37_chromosomes) > 1 or len(v38_chromosomes) > 1:
            error_msg = (
                f"Rs id infos for gene {gene} disagree on chromosome:\n"
                f"v37 chromosomes: {v37_chromosomes}\n"
                f"v38 chromosomes: {v38_chromosomes}"
            )
            raise ValueError(error_msg)

    @staticmethod
    def __assert_variants_in_haplotypes_compatible_with_rs_id_infos(
        haplotypes: FrozenSet[Haplotype], rs_id_infos: FrozenSet[RsIdInfo]
    ) -> None:
        variants = {variant for haplotype in haplotypes for variant in haplotype.variants}
        for variant in variants:
            matching_rs_id_infos = [rs_id_info for rs_id_info in rs_id_infos if rs_id_info.rs_id == variant.rs_id]
            assert len(matching_rs_id_infos) == 1, (
                f"Unexpected number of rs id infos match rs id with variant from haplotype:\n"
                f"variant={variant}, matches={matching_rs_id_infos}"
            )
            if variant.variant_allele == matching_rs_id_infos[0].reference_site_v38.allele:
                error_msg = (
                    f"Allele of variant matches reference allele from rs id info:\n"
                    f"variant={variant}, rs_id_info={matching_rs_id_infos[0]}"
                )
                raise ValueError(error_msg)

    @staticmethod
    def __assert_rs_ids_all_different(rs_id_infos: FrozenSet[RsIdInfo]) -> None:
        rs_ids = [info.rs_id for info in rs_id_infos]
        if len(rs_ids) != len(set(rs_ids)):
            error_msg = (
                f"Not all rs ids are different: rs_ids={sorted(rs_ids)}"
            )
            raise ValueError(error_msg)


def assert_no_overlap_gene_names(gene_infos: Collection[GeneInfo], source_name: str) -> None:
    if gene_names_overlap(gene_infos):
        gene_name_to_multiple_infos = get_gene_name_to_multiple_infos(gene_infos)
        error_msg = (
            f"The {source_name} contains gene summaries with the same gene names but different contents. "
            f"Duplicates: {gene_name_to_multiple_infos}"
        )
        raise ValueError(error_msg)


def gene_names_overlap(gene_infos: Collection[GeneInfo]) -> bool:
    return len({info.gene for info in gene_infos}) != len(gene_infos)


def get_gene_name_to_multiple_infos(gene_infos: Collection[GeneInfo]) -> Dict[str, List[GeneInfo]]:
    return get_key_to_multiple_values([(info.gene, info) for info in gene_infos])