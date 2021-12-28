from typing import FrozenSet, Optional, Set, Dict

from panel.drug_summary import DrugSummary
from panel.gene_transcript_summary import GeneTranscriptSummary
from panel.haplotype import Haplotype, GeneHaplotypePanel
from panel.rs_id_info import RsIdInfo
from panel.variant import Variant


class GenePanel(object):
    """This object is meant to be immutable"""

    def __init__(
        self,
        gene: str,
        wild_type_haplotype_name: str,
        transcript_id: Optional[str],
        other_haplotypes: FrozenSet[Haplotype],
        rs_id_infos: FrozenSet[RsIdInfo],
        drug_summaries: FrozenSet[DrugSummary],
    ) -> None:
        self.__assert_info_exists_for_all_rs_ids_in_haplotypes(other_haplotypes, rs_id_infos)
        self.__assert_variants_in_haplotypes_compatible_with_rs_id_infos(other_haplotypes, rs_id_infos)

        drug_name_to_drug_summary: Dict[str, DrugSummary] = {}
        for drug_summary in drug_summaries:
            if drug_summary.name in drug_name_to_drug_summary.keys():
                error_msg = f"The gene '{gene}' has multiple drugs with the name '{drug_summary.name}'."
                raise ValueError(error_msg)
            drug_name_to_drug_summary[drug_summary.name] = drug_summary

        self.__drug_name_to_drug_summary = drug_name_to_drug_summary
        self.__gene_haplotype_panel = GeneHaplotypePanel(gene, wild_type_haplotype_name, other_haplotypes)
        self.__gene_transcript_summary = GeneTranscriptSummary(gene, transcript_id, rs_id_infos)

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, GenePanel)
                and self.__drug_name_to_drug_summary == other.__drug_name_to_drug_summary
                and self.__gene_haplotype_panel == other.__gene_haplotype_panel
                and self.__gene_transcript_summary == other.__gene_transcript_summary
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.__get_drug_summaries(),
                self.__gene_haplotype_panel,
                self.__gene_transcript_summary,
            )
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"GeneInfo("
            f"gene={self.__gene_transcript_summary.gene!r}, "
            f"wild_type_haplotype_name={self.__gene_haplotype_panel.wild_type_haplotype_name!r}, "
            f"transcript_id={self.__gene_transcript_summary.transcript_id!r}, "
            f"haplotypes={self.__gene_haplotype_panel.get_haplotypes() !r}, "
            f"rs_id_infos={self.__gene_transcript_summary.rs_id_infos!r}, "
            f"drug_summaries={self.__get_drug_summaries()!r}, "
            f")"
        )

    @property
    def gene(self) -> str:
        return self.__gene_transcript_summary.gene

    @property
    def wild_type_haplotype_name(self) -> str:
        return self.__gene_haplotype_panel.wild_type_haplotype_name

    @property
    def transcript_id(self) -> Optional[str]:
        return self.__gene_transcript_summary.transcript_id

    @property
    def rs_id_infos(self) -> FrozenSet[RsIdInfo]:
        return self.__gene_transcript_summary.rs_id_infos

    def get_non_wild_type_haplotype_names(self) -> Set[str]:
        return self.__gene_haplotype_panel.get_non_wild_type_haplotype_names()

    def get_variants(self, haplotype_name: str) -> Set[Variant]:
        return self.__gene_haplotype_panel.get_variants(haplotype_name)

    def get_haplotype_function(self, haplotype_name: str) -> str:
        return self.__gene_haplotype_panel.get_haplotype_function(haplotype_name)

    def get_drug_names(self) -> Set[str]:
        return set(self.__drug_name_to_drug_summary.keys())

    def get_prescription_url(self, drug_name: str) -> str:
        return self.__drug_name_to_drug_summary[drug_name].url_prescription_info

    def get_rs_ids(self) -> Set[str]:
        return self.__gene_transcript_summary.get_rs_ids()

    def __get_drug_summaries(self) -> FrozenSet[DrugSummary]:
        return frozenset(self.__drug_name_to_drug_summary.values())

    @staticmethod
    def __assert_info_exists_for_all_rs_ids_in_haplotypes(
        non_wild_type_haplotypes: FrozenSet[Haplotype], rs_id_infos: FrozenSet[RsIdInfo]
    ) -> None:
        rs_ids_in_haplotypes = {
            variant.rs_id for haplotype in non_wild_type_haplotypes for variant in haplotype.variants
        }
        rs_ids_with_info = {info.rs_id for info in rs_id_infos}
        if not rs_ids_in_haplotypes.issubset(rs_ids_with_info):
            rs_ids_without_info = rs_ids_in_haplotypes.difference(rs_ids_with_info)
            error_msg = f"No info available for some of the rs ids in the haplotypes. Rs ids: {rs_ids_without_info}"
            raise ValueError(error_msg)

    @staticmethod
    def __assert_variants_in_haplotypes_compatible_with_rs_id_infos(
            non_wild_type_haplotypes: FrozenSet[Haplotype], rs_id_infos: FrozenSet[RsIdInfo]
    ) -> None:
        variants = {variant for haplotype in non_wild_type_haplotypes for variant in haplotype.variants}
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
