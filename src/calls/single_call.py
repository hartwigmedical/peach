from typing import NamedTuple, Tuple, Optional, FrozenSet

from calls.dual_call import DualCall
from util.filter import VcfCallFilter
from util.reference_assembly import ReferenceAssembly
from util.reference_site import ReferenceSite


class AnnotatedSingleCall(NamedTuple):
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


class VcfCall(NamedTuple):
    """
    Call as it was extracted from the input VCF.
    Has data and annotation for only a single reference assembly version. So only v37 or v38.
    """
    reference_site: ReferenceSite
    alleles: Tuple[str, str]  # The order is (ref, alt) when there is one of each
    rs_ids: Tuple[str, ...]
    variant_annotation: Optional[str]  # Is None if unknown
    filter: VcfCallFilter

    def matches(self, call: DualCall, reference_assembly: ReferenceAssembly) -> bool:
        if reference_assembly == ReferenceAssembly.V37:
            result = (
                    self.reference_site == call.reference_site_v37
                    and set(self.alleles).issubset({call.reference_site_v37.allele, call.alt_allele_v37})
            )
        elif reference_assembly == ReferenceAssembly.V38:
            result = (
                    self.reference_site == call.reference_site_v38
                    and set(self.alleles).issubset({call.reference_site_v38.allele, call.alt_allele_v38})
            )
        else:
            raise ValueError(f"Unexpected reference assembly: {reference_assembly}")
        return result

class VcfCallData(NamedTuple):
    calls: FrozenSet[VcfCall]
    reference_assembly: ReferenceAssembly

    def __repr__(self) -> str:  # pragma: no cover
        calls_string = "frozenset(" + ", ".join(sorted([repr(call) for call in self.calls])) + ")"
        return f"VcfCallData({calls_string}, {self.reference_assembly!r})"
