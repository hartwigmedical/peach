import unittest
from typing import Dict, FrozenSet

from base.gene_coordinate import GeneCoordinate
from base.reference_site import ReferenceSite
from config.annotation import Annotation
from config.drug_info import DrugInfo
from config.gene_info import GeneInfo
from config.haplotype import Haplotype
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant
from json_parser import JsonParser
from main import load_panel
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
            "panelName": "fake_panel",
            "panelVersion": "0.3",
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
                        },
                        {
                            "positionV37": "98348885",
                            "positionV38": "97883329",
                            "rsid": "rs1801265",
                            "referenceAlleleV37": "G",
                            "referenceAlleleV38": "A"
                        }
                    ],
                    "drugs": [
                        {
                            "name": "5-Fluorouracil",
                            "urlPrescriptionInfo": "https://www.source_url.org/5-Fluorouracil"
                        },
                        {
                            "name": "Capecitabine",
                            "urlPrescriptionInfo": "https://www.source_url.org/Capecitabine"
                        }
                    ],
                    "refSeqDifferenceAnnotations": [
                        {
                            "rsid": "rs72549303",
                            "annotationV37": "6744CA>GA",
                            "annotationV38": "6744GA>CA"
                        },
                        {
                            "rsid": "rs1801265",
                            "annotationV37": "85C>T",
                            "annotationV38": "85T>C"
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
                            "urlPrescriptionInfo": "https://www.source_url.org/Aspirin"
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
                            "urlPrescriptionInfo": "https://www.source_url.org/Aspirin"
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
        panel = JsonParser.get_panel(panel_json)
        panel_expected = self.__get_expected_panel()
        self.assertEqual(panel_expected, panel)

    def test_panel_name_from_json(self) -> None:
        """Test extraction and formatting of panel name with version"""
        panel_path = get_test_resource("test_panel.json")
        panel = load_panel(str(panel_path))
        expected_panel_id = "fake_panel_v0.3"
        self.assertEqual(expected_panel_id, panel.get_id())

    def __get_expected_panel(self) -> Panel:
        dpyd_two_a_variant = Variant("rs3918290", "T")
        dpyd_two_b_variant = Variant("rs1801159", "C")
        dpyd_three_variant = Variant("rs72549303", "TG")
        fake_variant = Variant("rs1212125", "C")
        fake2_variant = Variant("rs1212127", "C")
        dpyd_haplotypes_expected = frozenset(
            {
                Haplotype("*2A", "No Function", frozenset({dpyd_two_a_variant})),
                Haplotype("*2B", "No Function", frozenset({dpyd_two_a_variant, dpyd_two_b_variant})),
                Haplotype("*3", "Normal Function", frozenset({dpyd_three_variant})),
            }
        )
        dpyd_rs_id_infos_expected = frozenset(
            {
                RsIdInfo(
                    "rs3918290",
                    ReferenceSite(GeneCoordinate("1", 97915614), "C"),
                    ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ),
                RsIdInfo(
                    "rs72549309",
                    ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"),
                    ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ),
                RsIdInfo(
                    "rs1801159",
                    ReferenceSite(GeneCoordinate("1", 97981395), "T"),
                    ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ),
                RsIdInfo(
                    "rs72549303",
                    ReferenceSite(GeneCoordinate("1", 97915621), "TG"),
                    ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ),
                RsIdInfo(
                    "rs1801265",
                    ReferenceSite(GeneCoordinate("1", 98348885), "G"),
                    ReferenceSite(GeneCoordinate("chr1", 97883329), "A"),
                ),
            }
        )
        dpyd_drugs_expected = frozenset(
            {
                DrugInfo("5-Fluorouracil", "https://www.source_url.org/5-Fluorouracil"),
                DrugInfo("Capecitabine", "https://www.source_url.org/Capecitabine"),
            }
        )
        dpyd_rs_id_to_difference_annotations = {
            "rs72549303": Annotation("6744CA>GA", "6744GA>CA"),
            "rs1801265": Annotation("85C>T", "85T>C"),
        }
        dpyd_gene_info_expected = GeneInfo(
            "DPYD",
            "*1",
            "ENST00000370192",
            dpyd_haplotypes_expected,
            dpyd_rs_id_infos_expected,
            dpyd_drugs_expected,
            dpyd_rs_id_to_difference_annotations,
        )
        fake_haplotypes_expected = frozenset(
            {
                Haplotype("*4A", "Reduced Function", frozenset({fake_variant})),
            }
        )
        fake_rs_id_infos_expected = frozenset(
            {
                RsIdInfo(
                    "rs1212125",
                    ReferenceSite(GeneCoordinate("5", 97915617), "T"),
                    ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ),
            }
        )
        fake_drugs_expected = frozenset(
            {
                DrugInfo("Aspirin", "https://www.source_url.org/Aspirin"),
            }
        )
        fake_rs_id_to_difference_annotations: Dict[str, Annotation] = {}
        fake_gene_info_expected = GeneInfo(
            "FAKE",
            "*1",
            None,
            fake_haplotypes_expected,
            fake_rs_id_infos_expected,
            fake_drugs_expected,
            fake_rs_id_to_difference_annotations,
        )
        fake2_haplotypes_expected = frozenset(
            {
                Haplotype("*4A", "Reduced Function", frozenset({fake2_variant})),
            }
        )
        fake2_rs_id_infos_expected = frozenset(
            {
                RsIdInfo(
                    "rs1212127",
                    ReferenceSite(GeneCoordinate("16", 97915617), "C"),
                    ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ),
            }
        )
        fake2_drugs_expected = frozenset(
            {
                DrugInfo("Aspirin", "https://www.source_url.org/Aspirin"),
            }
        )
        fake2_rs_id_to_difference_annotations = {"rs1212127": Annotation("1324C>T", "1324T>C")}
        fake2_gene_info_expected = GeneInfo(
            "FAKE2",
            "*1",
            None,
            fake2_haplotypes_expected,
            fake2_rs_id_infos_expected,
            fake2_drugs_expected,
            fake2_rs_id_to_difference_annotations,
        )
        gene_infos_expected = frozenset(
            {
                dpyd_gene_info_expected,
                fake_gene_info_expected,
                fake2_gene_info_expected,
            }
        )
        name_expected = "fake_panel"
        version_expected = "0.3"
        panel_expected = Panel(name_expected, version_expected, gene_infos_expected)
        return panel_expected


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
