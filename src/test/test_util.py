import unittest
from typing import Set

from util.gene_coordinate import GeneCoordinate
from util.reference_site import ReferenceSite


class TestUtil(unittest.TestCase):
    def test_get_covered_coordinates_empty(self) -> None:
        reference = ReferenceSite(GeneCoordinate("X", 17), "")
        result = reference.get_covered_coordinates()

        result_expected: Set[GeneCoordinate] = set()
        self.assertEqual(result_expected, result)

    def test_get_covered_coordinates_single(self) -> None:
        reference = ReferenceSite(GeneCoordinate("X", 17), "A")
        result = reference.get_covered_coordinates()

        result_expected = {GeneCoordinate("X", 17)}
        self.assertEqual(result_expected, result)

    def test_get_covered_coordinates_multiple(self) -> None:
        reference = ReferenceSite(GeneCoordinate("X", 17), "AGTA")
        result = reference.get_covered_coordinates()

        result_expected = {
            GeneCoordinate("X", 17),
            GeneCoordinate("X", 18),
            GeneCoordinate("X", 19),
            GeneCoordinate("X", 20),
        }
        self.assertEqual(result_expected, result)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
