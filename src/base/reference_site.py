from typing import NamedTuple, Set

from base.gene_coordinate import GeneCoordinate


class ReferenceSite(NamedTuple):
    start_coordinate: GeneCoordinate
    allele: str

    def get_covered_coordinates(self) -> Set[GeneCoordinate]:
        covered_coordinates = {
            GeneCoordinate(self.start_coordinate.chromosome, self.start_coordinate.position + i)
            for i in range(len(self.allele))
        }
        return covered_coordinates
