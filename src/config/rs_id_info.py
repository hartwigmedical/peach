from typing import NamedTuple, Set

from base.gene_coordinate import GeneCoordinate
from base.json_alias import Json
from base.reference_assembly import ReferenceAssembly
from base.util import get_covered_coordinates


class RsIdInfo(NamedTuple):
    rs_id: str
    reference_allele_v37: str
    reference_allele_v38: str
    start_coordinate_v37: GeneCoordinate
    start_coordinate_v38: GeneCoordinate

    @classmethod
    def from_json(cls, data: Json, chromosome_v37: str, chromosome_v38: str) -> "RsIdInfo":
        rs_id = str(data['rsid'])
        reference_allele_v37 = str(data['referenceAlleleV37'])
        reference_allele_v38 = str(data['referenceAlleleV38'])
        start_coordinate_v37 = GeneCoordinate(chromosome_v37, int(data['positionV37']))
        start_coordinate_v38 = GeneCoordinate(chromosome_v38, int(data['positionV38']))
        info = RsIdInfo(
            rs_id,
            reference_allele_v37,
            reference_allele_v38,
            start_coordinate_v37,
            start_coordinate_v38,
        )
        return info

    def is_compatible(self, other: "RsIdInfo") -> bool:
        if self.rs_id == other.rs_id:
            return self == other
        else:
            v37_coordinates_overlap = self.get_relevant_coordinates(ReferenceAssembly.V37).intersection(
                other.get_relevant_coordinates(ReferenceAssembly.V37))
            v38_coordinates_overlap = self.get_relevant_coordinates(ReferenceAssembly.V38).intersection(
                other.get_relevant_coordinates(ReferenceAssembly.V38))
            return not v37_coordinates_overlap and not v38_coordinates_overlap

    def get_start_coordinate(self, reference_assembly: ReferenceAssembly) -> GeneCoordinate:
        if reference_assembly == ReferenceAssembly.V37:
            return self.start_coordinate_v37
        elif reference_assembly == ReferenceAssembly.V38:
            return self.start_coordinate_v38
        else:
            error_msg = "Unrecognized reference assembly version"
            raise NotImplementedError(error_msg)

    def get_reference_allele(self, reference_assembly: ReferenceAssembly) -> str:
        if reference_assembly == ReferenceAssembly.V37:
            return self.reference_allele_v37
        elif reference_assembly == ReferenceAssembly.V38:
            return self.reference_allele_v38
        else:
            error_msg = "Unrecognized reference assembly version"
            raise NotImplementedError(error_msg)

    def get_relevant_coordinates(self, reference_assembly: ReferenceAssembly) -> Set[GeneCoordinate]:
        start_coordinate = self.get_start_coordinate(reference_assembly)
        reference_allele = self.get_reference_allele(reference_assembly)
        return get_covered_coordinates(start_coordinate, reference_allele)
