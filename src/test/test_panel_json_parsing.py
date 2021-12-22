import unittest

from panel.panel import Panel
from json_parser import JsonParser
from main import load_panel
from test.util_for_test import get_wide_example_panel
from test_resources.test_resource import get_test_resource


class TestPanelJsonParsing(unittest.TestCase):
    def test_load_panel(self) -> None:
        """Load panel from json file"""
        panel_path = get_test_resource("test_panel.json")
        panel = load_panel(str(panel_path))

        panel_expected = self.__get_expected_panel()

        self.assertEqual(panel_expected, panel)

    def test_json_parser(self) -> None:
        """Load panel from json dictionary"""
        panel_json = {
            "panelName": "WideTestPanel",
            "panelVersion": "1.0",
            "genes": [
                {
                    "gene": "DPYD",
                    "chromosomeV37": "1",
                    "chromosomeV38": "chr1",
                    "wildTypeHaplotype": "*1",
                    "canonicalTranscript": "ENST00000370192",
                    "haplotypes": [
                        {
                            "haplotypeName": "*2A",
                            "function": "No Function",
                            "haplotypeVariants": [
                                {
                                    "rsid": "rs3918290",
                                    "altAlleleV38": "T"
                                }
                            ]
                        },
                        {
                            "haplotypeName": "*2B",
                            "function": "No Function",
                            "haplotypeVariants": [
                                {
                                    "rsid": "rs3918290",
                                    "altAlleleV38": "T"
                                },
                                {
                                    "rsid": "rs1801159",
                                    "altAlleleV38": "C"
                                }
                            ]
                        },
                        {
                            "haplotypeName": "*3",
                            "function": "Normal Function",
                            "haplotypeVariants": [
                                {
                                    "rsid": "rs72549303",
                                    "altAlleleV38": "TG"
                                }
                            ]
                        }
                    ],
                    "variants": [
                        {
                            "positionV37": "97915614",
                            "positionV38": "97450058",
                            "rsid": "rs3918290",
                            "referenceAlleleV37": "C",
                            "referenceAlleleV38": "C"
                        },
                        {
                            "positionV37": "98205966",
                            "positionV38": "97740410",
                            "rsid": "rs72549309",
                            "referenceAlleleV37": "GATGA",
                            "referenceAlleleV38": "GATGA"
                        },
                        {
                            "positionV37": "97981395",
                            "positionV38": "97515839",
                            "rsid": "rs1801159",
                            "referenceAlleleV37": "T",
                            "referenceAlleleV38": "T"
                        },
                        {
                            "positionV37": "97915621",
                            "positionV38": "97450065",
                            "rsid": "rs72549303",
                            "referenceAlleleV37": "TG",
                            "referenceAlleleV38": "TC"
                        }
                    ],
                    "drugs": [
                        {
                            "name": "5-Fluorouracil",
                            "urlPrescriptionInfo": "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"
                        },
                        {
                            "name": "Capecitabine",
                            "urlPrescriptionInfo": "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"
                        }
                    ],
                    "refSeqDifferenceAnnotations": [
                        {
                            "rsid": "rs72549303",
                            "annotationV37": "6744CA>GA",
                            "annotationV38": "6744GA>CA"
                        }
                    ]
                },
                {
                    "gene": "FAKE",
                    "chromosomeV37": "5",
                    "chromosomeV38": "chr5",
                    "wildTypeHaplotype": "*1",
                    "haplotypes": [
                        {
                            "haplotypeName": "*4A",
                            "function": "Reduced Function",
                            "haplotypeVariants": [
                                {
                                    "rsid": "rs1212125",
                                    "altAlleleV38": "C"
                                }
                            ]
                        }
                    ],
                    "variants": [
                        {
                            "positionV37": "97915617",
                            "positionV38": "97450060",
                            "rsid": "rs1212125",
                            "referenceAlleleV37": "T",
                            "referenceAlleleV38": "T"
                        }
                    ],
                    "drugs": [
                        {
                            "name": "Aspirin",
                            "urlPrescriptionInfo": "https://www.pharmgkb.org/some_other_url"
                        }
                    ],
                    "refSeqDifferenceAnnotations": [
                    ]
                },
                {
                    "gene": "FAKE2",
                    "chromosomeV37": "16",
                    "chromosomeV38": "chr16",
                    "wildTypeHaplotype": "*1",
                    "haplotypes": [
                        {
                            "haplotypeName": "*4A",
                            "function": "Reduced Function",
                            "haplotypeVariants": [
                                {
                                    "rsid": "rs1212127",
                                    "altAlleleV38": "C"
                                }
                            ]
                        }
                    ],
                    "variants": [
                        {
                            "positionV37": "97915617",
                            "positionV38": "97450060",
                            "rsid": "rs1212127",
                            "referenceAlleleV37": "C",
                            "referenceAlleleV38": "T"
                        }
                    ],
                    "drugs": [
                        {
                            "name": "Aspirin",
                            "urlPrescriptionInfo": "https://www.pharmgkb.org/some_other_url"
                        }
                    ],
                    "refSeqDifferenceAnnotations": [
                        {
                            "rsid": "rs1212127",
                            "annotationV37": "1324C>T",
                            "annotationV38": "1324T>C"
                        }
                    ]
                }
            ]
        }
        panel = JsonParser().get_panel(panel_json)
        panel_expected = self.__get_expected_panel()

        self.assertEqual(panel_expected, panel)

    def test_panel_name_from_json(self) -> None:
        """Test extraction and formatting of panel name with version"""
        panel_path = get_test_resource("test_panel.json")
        panel = load_panel(str(panel_path))
        expected_panel_id = "WideTestPanel_v1.0"
        self.assertEqual(expected_panel_id, panel.get_id())

    def __get_expected_panel(self) -> Panel:
        return get_wide_example_panel()


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
