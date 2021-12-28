from typing import Optional, NamedTuple, Tuple, FrozenSet

from base.filter import FullCallFilter
from base.reference_assembly import ReferenceAssembly
from base.reference_site import ReferenceSite


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

    def get_reference_site(self, reference_assembly: ReferenceAssembly) -> Optional[ReferenceSite]:
        if reference_assembly == ReferenceAssembly.V37:
            return self.reference_site_v37
        elif reference_assembly == ReferenceAssembly.V38:
            return self.reference_site_v38
        else:
            raise NotImplementedError(f"Unrecognized reference assembly: {reference_assembly}")


class FullCallData(NamedTuple):
    calls: FrozenSet[FullCall]

    def __repr__(self) -> str:  # pragma: no cover
        calls_string = ", ".join(sorted([repr(call) for call in self.calls]))
        return f"FullCallData(frozenset({calls_string}))"
