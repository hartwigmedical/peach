from typing import Collection, Dict, List, FrozenSet, Set

from base.util import get_key_to_multiple_values
from panel.variant import Variant


class Haplotype(object):
    def __init__(self, name: str, function: str, variants: FrozenSet[Variant]) -> None:
        if not variants:
            raise ValueError("Haplotype without variants is not allowed")

        seen_rs_ids: Set[str] = set()
        for variant in variants:
            if variant.rs_id in seen_rs_ids:
                error_msg = (
                    f"The haplotype '{name}' contains multiple variants with the same rs_id '{variant.rs_id}'."
                )
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


def assert_no_overlap_haplotype_names(haplotypes: Collection[Haplotype], source_name: str) -> None:
    if names_of_haplotypes_overlap(haplotypes):
        name_to_multiple_haplotypes = get_name_to_multiple_haplotypes(haplotypes)
        error_msg = (
            f"The {source_name} contains haplotypes with the same name but different summaries. "
            f"Duplicates: {name_to_multiple_haplotypes}"
        )
        raise ValueError(error_msg)


def names_of_haplotypes_overlap(haplotypes: Collection[Haplotype]) -> bool:
    return len({haplotype.name for haplotype in haplotypes}) != len(haplotypes)


def get_name_to_multiple_haplotypes(haplotypes: Collection[Haplotype]) -> Dict[str, List[Haplotype]]:
    return get_key_to_multiple_values([(haplotype.name, haplotype) for haplotype in haplotypes])


def assert_no_overlap_haplotype_variant_combinations(haplotypes: Collection[Haplotype], source_name: str) -> None:
    if variant_combinations_of_haplotypes_overlap(haplotypes):
        variant_combination_to_multiple_haplotypes = get_variant_combination_to_multiple_haplotypes(haplotypes)
        error_msg = (
            f"The {source_name} contains haplotypes with the same variant combination but different names. "
            f"Duplicates: {variant_combination_to_multiple_haplotypes}"
        )
        raise ValueError(error_msg)


def variant_combinations_of_haplotypes_overlap(haplotypes: Collection[Haplotype]) -> bool:
    return len({haplotype.variants for haplotype in haplotypes}) != len(haplotypes)


def get_variant_combination_to_multiple_haplotypes(
    haplotypes: Collection[Haplotype],
) -> Dict[FrozenSet[Variant], List[Haplotype]]:
    return get_key_to_multiple_values([(haplotype.variants, haplotype) for haplotype in haplotypes])
