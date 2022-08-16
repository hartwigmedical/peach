import logging
from typing import Collection, Dict, FrozenSet, Set

from util.constants import NORMAL_FUNCTION_STRING
from panel.variant import Variant


class Haplotype(object):
    def __init__(self, name: str, function: str, variants: FrozenSet[Variant]) -> None:
        if not variants:
            raise ValueError("Haplotype without variants is not allowed")

        seen_rs_ids: Set[str] = set()
        for variant in variants:
            if variant.rs_id in seen_rs_ids:
                error_msg = f"The haplotype '{name}' contains multiple variants with the same rs_id '{variant.rs_id}'."
                raise ValueError(error_msg)
            seen_rs_ids.add(variant.rs_id)

        self.__name = name
        self.__function = function
        self.__variants = variants

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Haplotype)
            and self.__name == other.__name
            and self.__function == other.__function
            and self.__variants == other.__variants
        )

    def __hash__(self) -> int:
        return hash((self.__name, self.__function, self.__variants))

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"Haplotype("
            f"name={self.__name!r}, "
            f"function={self.__function!r}, "
            f"variants={self.__variants!r}, "
            f")"
        )

    @property
    def name(self) -> str:
        return self.__name

    @property
    def function(self) -> str:
        return self.__function

    @property
    def variants(self) -> FrozenSet[Variant]:
        return self.__variants


class GeneHaplotypePanel(object):
    def __init__(self, gene: str, wild_type_haplotype_name: str, other_haplotypes: FrozenSet[Haplotype]) -> None:
        self.__assert_no_overlap_haplotype_variant_combinations(other_haplotypes, gene)

        if not other_haplotypes:
            logging.warning(f"No alternate haplotypes configured for gene {gene}\n")

        haplotype_name_to_haplotype: Dict[str, Haplotype] = {}
        for haplotype in other_haplotypes:
            if haplotype.name in haplotype_name_to_haplotype.keys():
                error_msg = f"The gene '{gene}' has multiple haplotypes with the name '{haplotype.name}'."
                raise ValueError(error_msg)
            haplotype_name_to_haplotype[haplotype.name] = haplotype

        # public through @property decorator
        self.__gene = gene
        self.__wild_type_haplotype_name = wild_type_haplotype_name

        # truly private
        self.__haplotype_name_to_haplotype = haplotype_name_to_haplotype

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, GeneHaplotypePanel)
                and self.__gene == other.__gene
                and self.__wild_type_haplotype_name == other.__wild_type_haplotype_name
                and self.__haplotype_name_to_haplotype == other.__haplotype_name_to_haplotype
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.__gene,
                self.__wild_type_haplotype_name,
                self.get_haplotypes(),
            )
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"GeneHaplotypePanel("
            f"gene={self.__gene!r}, "
            f"wild_type_haplotype_name={self.__wild_type_haplotype_name!r}, "
            f"haplotypes={self.get_haplotypes()!r}, "
            f")"
        )

    @property
    def gene(self) -> str:
        return self.__gene

    @property
    def wild_type_haplotype_name(self) -> str:
        return self.__wild_type_haplotype_name

    def get_non_wild_type_haplotype_names(self) -> Set[str]:
        return set(self.__haplotype_name_to_haplotype.keys())

    def get_variants(self, haplotype_name: str) -> Set[Variant]:
        if haplotype_name == self.__wild_type_haplotype_name:
            return set()
        else:
            return set(self.__haplotype_name_to_haplotype[haplotype_name].variants)

    def get_haplotype_function(self, haplotype_name: str) -> str:
        if haplotype_name == self.__wild_type_haplotype_name:
            return NORMAL_FUNCTION_STRING
        else:
            return self.__haplotype_name_to_haplotype[haplotype_name].function

    def get_haplotypes(self) -> FrozenSet[Haplotype]:
        return frozenset(self.__haplotype_name_to_haplotype.values())

    @classmethod
    def __assert_no_overlap_haplotype_variant_combinations(
            cls, haplotypes: Collection[Haplotype], gene: str
    ) -> None:
        variant_combinations_of_haplotypes_overlap = (
            len({haplotype.variants for haplotype in haplotypes}) != len(haplotypes)
        )
        if variant_combinations_of_haplotypes_overlap:
            error_msg = (
                f"The gene '{gene}' has haplotypes with the same variant combination but different names."
            )
            raise ValueError(error_msg)
