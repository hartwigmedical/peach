from typing import NamedTuple, Set

from base.gene_coordinate import GeneCoordinate
from base.json_alias import Json
from base.reference_site import ReferenceSite
from base.reference_assembly import ReferenceAssembly


class RsIdInfo(NamedTuple):
    rs_id: str
    reference_site_v37: ReferenceSite
    reference_site_v38: ReferenceSite

    @classmethod
    def from_json(cls, data: Json, chromosome_v37: str, chromosome_v38: str) -> "RsIdInfo":
        rs_id = str(data["rsid"])
        reference_allele_v37 = str(data["referenceAlleleV37"])
        reference_allele_v38 = str(data["referenceAlleleV38"])
        start_coordinate_v37 = GeneCoordinate(chromosome_v37, int(data["positionV37"]))
        start_coordinate_v38 = GeneCoordinate(chromosome_v38, int(data["positionV38"]))
        info = RsIdInfo(
            rs_id,
            ReferenceSite(start_coordinate_v37, reference_allele_v37),
            ReferenceSite(start_coordinate_v38, reference_allele_v38),
        )
        return info

    def is_compatible(self, other: "RsIdInfo") -> bool:
        # TODO: check layout after black
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
