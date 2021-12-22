from typing import NamedTuple, Tuple, Optional

from base.filter import VcfCallFilter
from base.reference_site import ReferenceSite


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
