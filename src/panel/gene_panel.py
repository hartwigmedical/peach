import itertools
from typing import FrozenSet, Optional, Set, Dict

from panel.drug_info import DrugInfo
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
        drugs: FrozenSet[DrugInfo],
    ) -> None:
        self.__assert_rs_ids_all_different(rs_id_infos)
        self.__assert_rs_id_infos_compatible(rs_id_infos)
        self.__assert_rs_id_infos_match_chromosome(rs_id_infos, gene)
        self.__assert_info_exists_for_all_rs_ids_in_haplotypes(other_haplotypes, rs_id_infos)
        self.__assert_variants_in_haplotypes_compatible_with_rs_id_infos(other_haplotypes, rs_id_infos)

        drug_name_to_drug_info: Dict[str, DrugInfo] = {}
        for drug_info in drugs:
            if drug_info.name in drug_name_to_drug_info.keys():
                error_msg = f"The gene '{gene}' has multiple drugs with the name '{drug_info.name}'."
                raise ValueError(error_msg)
            drug_name_to_drug_info[drug_info.name] = drug_info

        # public through @property decorator
        self.__gene = gene
        self.__transcript_id = transcript_id
        self.__rs_id_infos = rs_id_infos

        # truly private
        self.__drug_name_to_drug_info = drug_name_to_drug_info
        self.__gene_haplotype_panel = GeneHaplotypePanel(gene, wild_type_haplotype_name, other_haplotypes)

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, GenePanel)
            and self.__gene == other.__gene
            and self.__transcript_id == other.__transcript_id
            and self.__rs_id_infos == other.__rs_id_infos
            and self.__drug_name_to_drug_info == other.__drug_name_to_drug_info
            and self.__gene_haplotype_panel == other.__gene_haplotype_panel
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.__gene,
                self.__transcript_id,
                self.__rs_id_infos,
                self.__get_drug_infos(),
                self.__gene_haplotype_panel,
            )
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"GeneInfo("
            f"gene={self.__gene!r}, "
            f"wild_type_haplotype_name={self.__gene_haplotype_panel.wild_type_haplotype_name!r}, "
            f"transcript_id={self.__transcript_id!r}, "
            f"haplotypes={self.__gene_haplotype_panel.get_haplotypes() !r}, "
            f"rs_id_infos={self.__rs_id_infos!r}, "
            f"drugs={self.__get_drug_infos()!r}, "
            f")"
        )

    @property
    def gene(self) -> str:
        return self.__gene

    @property
    def wild_type_haplotype_name(self) -> str:
        return self.__gene_haplotype_panel.wild_type_haplotype_name

    @property
    def transcript_id(self) -> Optional[str]:
        return self.__transcript_id

    @property
    def rs_id_infos(self) -> FrozenSet[RsIdInfo]:
        return self.__rs_id_infos

    def get_non_wild_type_haplotype_names(self) -> Set[str]:
        return self.__gene_haplotype_panel.get_non_wild_type_haplotype_names()

    def get_variants(self, haplotype_name: str) -> Set[Variant]:
        return self.__gene_haplotype_panel.get_variants(haplotype_name)

    def get_haplotype_function(self, haplotype_name: str) -> str:
        return self.__gene_haplotype_panel.get_haplotype_function(haplotype_name)

    def get_drug_names(self) -> Set[str]:
        return set(self.__drug_name_to_drug_info.keys())

    def get_prescription_url(self, drug_name: str) -> str:
        return self.__drug_name_to_drug_info[drug_name].url_prescription_info

    def get_rs_ids(self) -> Set[str]:
        return {rs_id_info.rs_id for rs_id_info in self.__rs_id_infos}

    def __get_drug_infos(self) -> FrozenSet[DrugInfo]:
        return frozenset(self.__drug_name_to_drug_info.values())

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

    @staticmethod
    def __assert_rs_ids_all_different(rs_id_infos: FrozenSet[RsIdInfo]) -> None:
        rs_ids = [info.rs_id for info in rs_id_infos]
        if len(rs_ids) != len(set(rs_ids)):
            error_msg = (
                f"Not all rs ids are different: rs_ids={sorted(rs_ids)}"
            )
            raise ValueError(error_msg)
