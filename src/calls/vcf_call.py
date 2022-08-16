from typing import NamedTuple, Tuple, Optional, FrozenSet

from util.filter import VcfCallFilter
from util.reference_assembly import ReferenceAssembly
from util.reference_site import ReferenceSite


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
        return f"VcfCallData({calls_string}, {self.reference_assembly!r})"
