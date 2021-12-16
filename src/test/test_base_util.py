import unittest
from typing import List, Tuple, Dict, Set

from base.gene_coordinate import GeneCoordinate
from base.reference_site import ReferenceSite
from base.util import get_key_to_multiple_values, replace_file_extension_of_path, strip_prefix


class TestBaseUtil(unittest.TestCase):
    def test_get_key_to_multiple_values_empty(self) -> None:
        key_value_pairs: List[Tuple[str, int]] = []
        result = get_key_to_multiple_values(key_value_pairs)
        result_expected: Dict[str, List[int]] = {}
        self.assertEqual(result_expected, result)

    def test_get_key_to_multiple_values_non_empty(self) -> None:
        key_value_pairs = [("A", 1), ("A", 2), ("A", 3), ("B", 1), ("B", 5), ("C", 1), ("D", 9)]
        result = get_key_to_multiple_values(key_value_pairs)

        result_expected = {
            "A": [1, 2, 3],
            "B": [1, 5],
        }

        self.assertEqual(result_expected, result)

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

    def test_replace_file_extension_of_path_single_point(self) -> None:
        result = replace_file_extension_of_path("something.ext", "test")
        result_expected = "something.test"
        self.assertEqual(result_expected, result)

    def test_replace_file_extension_of_path_multiple_points(self) -> None:
        result = replace_file_extension_of_path("this.is.something.ext", "test")
        result_expected = "this.is.something.test"
        self.assertEqual(result_expected, result)

    def test_replace_file_extension_of_path_no_point(self) -> None:
        with self.assertRaises(AssertionError):
            replace_file_extension_of_path("something", "test")

    def test_strip_prefix(self) -> None:
        self.assertEqual("HiHello", strip_prefix("HiHiHello", "Hi"))

    def test_strip_prefix_wrong_prefix(self) -> None:
        with self.assertRaises(ValueError):
            strip_prefix("HiHiHello", "Hello")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
