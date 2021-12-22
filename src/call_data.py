from copy import deepcopy
from typing import NamedTuple, Tuple, Optional, FrozenSet, Dict

from base.filter import FullCallFilter, VcfCallFilter
from base.reference_site import ReferenceSite
from base.reference_assembly import ReferenceAssembly


class VcfCall(NamedTuple):
    """
    Call as it was extracted from the input VCF.
    Has data and annotation for only a single reference assembly version. So only v37 or v38.
    """
    reference_site: ReferenceSite
    alleles: Tuple[str, str]  # The order is (ref, alt) when there is one of each
    gene: Optional[str]  # Is None if unknown
    rs_ids: Tuple[str, ...]
    variant_annotation: Optional[str]  # Is None if unknown
    filter: VcfCallFilter


class VcfCallData(NamedTuple):
    calls: FrozenSet[VcfCall]
    reference_assembly: ReferenceAssembly

    def __repr__(self) -> str:  # pragma: no cover
        calls_string = "frozenset(" + ", ".join(sorted([repr(call) for call in self.calls])) + ")"
        return f"SimpleCallData({calls_string}, {self.reference_assembly!r})"


class SimpleCall(NamedTuple):
    """
    Call with data and annotation for only a single reference assembly version. So only v37 or v38.
    Annotation fixes etc. with data from the panel has already been done at this point.
    """
    reference_site: ReferenceSite
    alleles: Tuple[str, str]  # The order is (ref, alt) when there is one of each
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation: Optional[str]  # Is None if unknown
    filter: VcfCallFilter

    def is_pass(self) -> bool:
        if self.filter == VcfCallFilter.PASS:
            return True
        elif self.filter == VcfCallFilter.NO_CALL:
            return False
        else:
            raise NotImplementedError("Unrecognized filter value")


class AnnotatedAllele(object):
    def __init__(self, allele: str, reference_assembly_to_is_variant: Dict[ReferenceAssembly, bool]) -> None:
        self.__allele = allele
        self.__reference_assembly_to_is_variant = deepcopy(reference_assembly_to_is_variant)

    @classmethod
    def from_reference_sites(
        cls, allele: str, reference_assembly_to_reference_site: Dict[ReferenceAssembly, Optional[ReferenceSite]]
    ) -> "AnnotatedAllele":
        reference_assembly_to_is_variant: Dict[ReferenceAssembly, bool] = {}
        for reference_assembly, reference_site in reference_assembly_to_reference_site.items():
            if reference_site is not None:
                reference_assembly_to_is_variant[reference_assembly] = allele != reference_site.allele
        return AnnotatedAllele(allele, reference_assembly_to_is_variant)

    def __eq__(self, other: object) -> bool:  # pragma: no cover
        return (
            isinstance(other, AnnotatedAllele)
            and self.__allele == other.__allele
            and self.__reference_assembly_to_is_variant == other.__reference_assembly_to_is_variant
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"AnnotatedAllele("
            f"allele={self.__allele!r}, "
            f"reference_assembly_to_is_variant={self.__reference_assembly_to_is_variant!r}, "
            f")"
        )

    @property
    def allele(self) -> str:
        return self.__allele

    def is_annotated_vs(self, reference_assembly: ReferenceAssembly) -> bool:
        return reference_assembly in self.__reference_assembly_to_is_variant.keys()

    def is_variant_vs(self, reference_assembly: ReferenceAssembly) -> bool:
        try:
            return self.__reference_assembly_to_is_variant[reference_assembly]
        except KeyError:
            raise SyntaxError(f"Allele not annotated vs reference assembly: {reference_assembly}")


class FullCall(NamedTuple):
    # Call with both v37 and v38 data and annotation
    reference_site_v37: Optional[ReferenceSite]  # Is None if unknown
    reference_site_v38: Optional[ReferenceSite]  # Is None if unknown
    alleles: Tuple[str, str]
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation_v37: Optional[str]  # Is None if unknown
    filter_v37: FullCallFilter
    variant_annotation_v38: Optional[str]  # Is None if unknown
    filter_v38: FullCallFilter

    def get_annotated_alleles(self) -> Tuple[AnnotatedAllele, AnnotatedAllele]:
        return self.__annotate_allele(self.alleles[0]), self.__annotate_allele(self.alleles[1])

    def __annotate_allele(self, allele: str) -> AnnotatedAllele:
        reference_assembly_to_reference_allele = {
            ReferenceAssembly.V37: self.reference_site_v37,
            ReferenceAssembly.V38: self.reference_site_v38,
        }
        return AnnotatedAllele.from_reference_sites(allele, reference_assembly_to_reference_allele)


class FullCallData(NamedTuple):
    calls: FrozenSet[FullCall]

    def __repr__(self) -> str:  # pragma: no cover
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
        return f"HaplotypeCall({self.__haplotype_name!r}, {self.__count!r})"

    @property
    def haplotype_name(self) -> str:
        return self.__haplotype_name

    @property
    def count(self) -> int:
        return self.__count
