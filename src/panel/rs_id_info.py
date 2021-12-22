from typing import NamedTuple

from base.reference_site import ReferenceSite
from base.reference_assembly import ReferenceAssembly


class RsIdInfo(NamedTuple):
    rs_id: str
    reference_site_v37: ReferenceSite
    reference_site_v38: ReferenceSite

    def is_compatible(self, other: "RsIdInfo") -> bool:
        if self.rs_id == other.rs_id:
            return self == other
        else:
            v37_coordinates_overlap = (
                self.get_reference_site(ReferenceAssembly.V37)
                .get_covered_coordinates()
                .intersection(other.get_reference_site(ReferenceAssembly.V37).get_covered_coordinates())
            )
            v38_coordinates_overlap = (
                self.get_reference_site(ReferenceAssembly.V38)
                .get_covered_coordinates()
                .intersection(other.get_reference_site(ReferenceAssembly.V38).get_covered_coordinates())
            )
            return not v37_coordinates_overlap and not v38_coordinates_overlap

    def has_reference_sequence_difference(self) -> bool:
        return self.reference_site_v37.allele != self.reference_site_v38.allele

    def get_reference_site(self, reference_assembly: ReferenceAssembly) -> ReferenceSite:
        if reference_assembly == ReferenceAssembly.V37:
            return self.reference_site_v37
        elif reference_assembly == ReferenceAssembly.V38:
            return self.reference_site_v38
        else:
            error_msg = "Unrecognized reference assembly version"
            raise NotImplementedError(error_msg)
