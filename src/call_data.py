from typing import NamedTuple, Tuple, Optional, Set, FrozenSet

from base.filter import FullCallFilter, SimpleCallFilter
from base.gene_coordinate import GeneCoordinate
from base.reference_assembly import ReferenceAssembly
from base.util import get_covered_coordinates


class SimpleCall(NamedTuple):
    # Call with data and annotation for only a single reference assembly version. So only v37 or v38
    start_coordinate: GeneCoordinate
    reference_allele: str
    alleles: Tuple[str, str]  # The order is (ref, alt) when there is one of each
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation: str
    filter: SimpleCallFilter

    def get_relevant_coordinates(self) -> Set[GeneCoordinate]:
        return get_covered_coordinates(self.start_coordinate, self.reference_allele)

    def is_pass(self) -> bool:
        if self.filter == SimpleCallFilter.PASS:
            return True
        elif self.filter == SimpleCallFilter.NO_CALL:
            return False
        else:
            raise NotImplementedError("Unrecognized filer value")


class SimpleCallData(NamedTuple):
    calls: FrozenSet[SimpleCall]
    reference_assembly: ReferenceAssembly

    def __repr__(self) -> str:
        calls_string = "frozenset(" + ", ".join(sorted([repr(call) for call in self.calls])) + ")"
        return f"SimpleCallData({calls_string}, {self.reference_assembly!r})"


class AnnotatedAllele(object):
    def __init__(self, allele: str, is_variant_vs_v37: Optional[bool], is_variant_vs_v38: Optional[bool]) -> None:
        self.__allele = allele
        self.__is_variant_vs_v37 = is_variant_vs_v37  # is None if unknown
        self.__is_variant_vs_v38 = is_variant_vs_v38  # is None if unknown

    @classmethod
    def from_alleles(cls, allele: str, reference_allele_v37: Optional[str],
                     reference_allele_v38: Optional[str]) -> "AnnotatedAllele":
        is_v37_variant = (
            (allele != reference_allele_v37) if reference_allele_v37 is not None else None
        )
        is_v38_variant = (
            (allele != reference_allele_v38) if reference_allele_v38 is not None else None
        )
        return AnnotatedAllele(allele, is_v37_variant, is_v38_variant)

    def __eq__(self, other: object) -> bool:  # pragma: no cover
        return (
                isinstance(other, AnnotatedAllele)
                and self.__allele == other.__allele
                and self.__is_variant_vs_v37 == other.__is_variant_vs_v37
                and self.__is_variant_vs_v38 == other.__is_variant_vs_v38
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"AnnotatedAllele("
            f"allele={self.__allele!r}, "
            f"is_variant_vs_v37={self.__is_variant_vs_v37!r}, "
            f"is_variant_vs_v38={self.__is_variant_vs_v38!r}, "
            f")"
        )

    @property
    def allele(self) -> str:
        return self.__allele

    @property
    def is_variant_vs_v37(self) -> bool:
        if self.__is_variant_vs_v37 is None:
            raise ValueError("Cannot get is_variant_vs_v37 if it is None")
        return self.__is_variant_vs_v37

    def is_annotated_vs_v37(self) -> bool:
        return self.__is_variant_vs_v37 is not None

    @property
    def is_variant_vs_v38(self) -> bool:
        if self.__is_variant_vs_v38 is None:
            raise ValueError("Cannot get is_variant_vs_v38 if it is None")
        return self.__is_variant_vs_v38

    def is_annotated_vs_v38(self) -> bool:
        return self.__is_variant_vs_v38 is not None


class FullCall(NamedTuple):
    # Call with both v37 and v38 data and annotation
    start_coordinate_v37: Optional[GeneCoordinate]  # is None if unknown
    reference_allele_v37: Optional[str]  # is None if unknown
    start_coordinate_v38: Optional[GeneCoordinate]  # is None if unknown
    reference_allele_v38: Optional[str]  # is None if unknown
    alleles: Tuple[str, str]
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation_v37: str
    filter_v37: FullCallFilter
    variant_annotation_v38: str
    filter_v38: FullCallFilter

    def get_relevant_v37_coordinates(self) -> Optional[Set[GeneCoordinate]]:
        # Returns None if missing some information to determine these coordinates
        if self.start_coordinate_v37 is None or self.reference_allele_v37 is None:
            return None
        else:
            return get_covered_coordinates(self.start_coordinate_v37, self.reference_allele_v37)

    def get_relevant_v38_coordinates(self) -> Optional[Set[GeneCoordinate]]:
        # Returns None if missing some information to determine these coordinates
        if self.start_coordinate_v38 is None or self.reference_allele_v38 is None:
            return None
        else:
            return get_covered_coordinates(self.start_coordinate_v38, self.reference_allele_v38)

    def get_annotated_alleles(self) -> Tuple[AnnotatedAllele, AnnotatedAllele]:
        annotated_alleles = self.__annotate_allele(self.alleles[0]), self.__annotate_allele(self.alleles[1])
        self.__assert_alleles_in_expected_order(annotated_alleles)
        return annotated_alleles

    def __annotate_allele(self, allele: str) -> AnnotatedAllele:
        return AnnotatedAllele.from_alleles(allele, self.reference_allele_v37, self.reference_allele_v38)

    @classmethod
    def __assert_alleles_in_expected_order(cls, annotated_alleles: Tuple[AnnotatedAllele, AnnotatedAllele]) -> None:
        alleles_in_unexpected_order = (
                annotated_alleles[0].is_variant_vs_v37
                and not annotated_alleles[1].is_variant_vs_v37
        )
        if alleles_in_unexpected_order:
            error_msg = (f"Alleles are in unexpected order, alt before ref: "
                         f"alleles=({annotated_alleles[0]}, {annotated_alleles[1]})")
            raise ValueError(error_msg)


class FullCallData(NamedTuple):
    calls: FrozenSet[FullCall]

    def __repr__(self) -> str:
        calls_string = ", ".join(sorted([repr(call) for call in self.calls]))
        return f"FullCallData(frozenset({calls_string}))"


class HaplotypeCall(object):
    def __init__(self, haplotype_name: str, count: int) -> None:
        if not 1 <= count <= 2:
            error_msg = f"Illegal haplotype count {count} for haplotype {haplotype_name}"
            raise SyntaxError(error_msg)

        self.__haplotype_name = haplotype_name
        self.__count = count

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, HaplotypeCall)
                and self.__haplotype_name == other.__haplotype_name
                and self.__count == other.__count
        )

    def __hash__(self) -> int:
        return hash((self.__haplotype_name, self.__count))

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"HaplotypeCall("
            f"{self.__haplotype_name!r}, "
            f"{self.__count!r}, "
            f")"
        )

    @property
    def haplotype_name(self) -> str:
        return self.__haplotype_name

    @property
    def count(self) -> int:
        return self.__count
