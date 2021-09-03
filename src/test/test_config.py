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
from main import load_panel
from test_resources.test_resource import get_test_resource


class TestLoadConfig(unittest.TestCase):
    def test_load_panel(self) -> None:
        """Load panel from json"""
        panel_path = get_test_resource("test_panel.json")
        panel = load_panel(str(panel_path))

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

        self.assertEqual(panel_expected, panel)

    def test_panel_name_from_json(self) -> None:
        """Test extraction and formatting of panel name with version"""
        panel_path = get_test_resource("test_panel.json")
        panel = load_panel(str(panel_path))
        expected_panel_id = "fake_panel_v0.3"
        self.assertEqual(expected_panel_id, panel.get_id())

    def test_gene_info_with_overlapping_haplotype_names(self) -> None:
        """Error when haplotype name used multiple times for gene"""
        gene = "FAKE"
        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_diff_annotation: Dict[str, Annotation] = dict()

        variant1 = Variant("rs94982", "A")
        variant2 = Variant("rs394934", "T")

        variant1_rs_id_info = RsIdInfo(
            variant1.rs_id,
            ReferenceSite(GeneCoordinate(chromosome_v37, 4994545), "C"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 2993823), "C"),
        )
        variant2_rs_id_info = RsIdInfo(
            variant2.rs_id,
            ReferenceSite(GeneCoordinate(chromosome_v37, 3993842), "G"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 2949923), "G"),
        )
        rs_id_infos = frozenset(
            [
                variant1_rs_id_info,
                variant2_rs_id_info,
            ]
        )

        haplotypes1 = frozenset(
            [
                Haplotype("*2", "No Function", frozenset([variant1])),
                Haplotype("*4", "Partial Function", frozenset([variant1, variant2])),
            ]
        )
        haplotypes2 = frozenset(
            [
                Haplotype("*2", "No Function", frozenset([variant1])),
                Haplotype("*2", "Partial Function", frozenset([variant1, variant2])),
            ]
        )

        GeneInfo(gene, reference_haplotype_name, haplotypes1, rs_id_infos, drugs, rs_id_to_ref_seq_diff_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, reference_haplotype_name, haplotypes2, rs_id_infos, drugs, rs_id_to_ref_seq_diff_annotation)

    def test_gene_info_with_overlapping_haplotype_variants(self) -> None:
        """Error when different haplotypes have the exact same variant combination"""
        gene = "FAKE"
        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_diff_annotation: Dict[str, Annotation] = dict()

        variant1 = Variant("rs94982", "A")
        variant2 = Variant("rs394934", "T")
        variant3 = Variant("rs495825", "C")

        variant1_rs_id_info = RsIdInfo(
            variant1.rs_id,
            ReferenceSite(GeneCoordinate(chromosome_v37, 4994545), "C"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 2993823), "C"),
        )
        variant2_rs_id_info = RsIdInfo(
            variant2.rs_id,
            ReferenceSite(GeneCoordinate(chromosome_v37, 3993842), "G"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 2949923), "G"),
        )
        variant3_rs_id_info = RsIdInfo(
            variant3.rs_id,
            ReferenceSite(GeneCoordinate(chromosome_v37, 293923), "A"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 138812), "A"),
        )
        rs_id_infos = frozenset(
            [
                variant1_rs_id_info,
                variant2_rs_id_info,
                variant3_rs_id_info,
            ]
        )

        haplotypes1 = frozenset(
            [
                Haplotype("*2", "No Function", frozenset([variant1, variant2])),
                Haplotype("*4", "Partial Function", frozenset([variant1, variant3])),
            ]
        )
        haplotypes2 = frozenset(
            [
                Haplotype("*2", "No Function", frozenset([variant1, variant2])),
                Haplotype("*4", "Partial Function", frozenset([variant1, variant2])),
            ]
        )

        GeneInfo(gene, reference_haplotype_name, haplotypes1, rs_id_infos, drugs, rs_id_to_ref_seq_diff_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, reference_haplotype_name, haplotypes2, rs_id_infos, drugs, rs_id_to_ref_seq_diff_annotation)

    def test_gene_info_with_overlapping_drug_names(self) -> None:
        """Error when drug name used multiple times for gene"""
        gene = "FAKE"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        rs_id_infos: FrozenSet[RsIdInfo] = frozenset()
        rs_id_to_ref_seq_diff_annotation: Dict[str, Annotation] = dict()

        drug_info1 = DrugInfo("Paracetamol", "google.com")
        drug_info2 = DrugInfo("Paracetamol", "wikipedia.com")
        drugs = frozenset([drug_info1])
        overlapping_drugs = frozenset([drug_info1, drug_info2])

        GeneInfo(gene, reference_haplotype_name, haplotypes, rs_id_infos, drugs, rs_id_to_ref_seq_diff_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                haplotypes,
                rs_id_infos,
                overlapping_drugs,
                rs_id_to_ref_seq_diff_annotation,
            )

    def test_gene_info_with_overlapping_rs_id_infos(self) -> None:
        """Error when gene info has rs id infos for which the relevant coordinates overlap"""
        gene = "FAKE"
        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_diff_annotation: Dict[str, Annotation] = dict()

        rs_id_info1 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "A"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "A"),
        )
        rs_id_info2 = RsIdInfo(
            "rs294927",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499592), "AA"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399482), "AA"),
        )
        single_rs_id_info = frozenset([rs_id_info1])
        overlapping_rs_id_infos = frozenset([rs_id_info1, rs_id_info2])

        GeneInfo(gene, reference_haplotype_name, haplotypes, single_rs_id_info, drugs, rs_id_to_ref_seq_diff_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                haplotypes,
                overlapping_rs_id_infos,
                drugs,
                rs_id_to_ref_seq_diff_annotation,
            )

    def test_gene_info_with_rs_id_infos_for_different_chromosome(self) -> None:
        """Error when gene info has rs id infos for which the relevant coordinates have a different chromosome"""
        gene = "FAKE"
        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_diff_annotation: Dict[str, Annotation] = dict()

        other_chromosome_v37 = "1"
        other_chromosome_v38 = "1"

        rs_id_info1 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "A"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "A"),
        )
        rs_id_info2 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(other_chromosome_v37, 499593), "A"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "A"),
        )
        rs_id_info3 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "A"),
            ReferenceSite(GeneCoordinate(other_chromosome_v38, 399483), "A"),
        )
        rs_id_info4 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(other_chromosome_v37, 499593), "A"),
            ReferenceSite(GeneCoordinate(other_chromosome_v38, 399483), "A"),
        )

        GeneInfo(
            gene,
            reference_haplotype_name,
            haplotypes,
            frozenset([rs_id_info1]),
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        GeneInfo(
            gene,
            reference_haplotype_name,
            haplotypes,
            frozenset([rs_id_info2]),
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        GeneInfo(
            gene,
            reference_haplotype_name,
            haplotypes,
            frozenset([rs_id_info3]),
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        GeneInfo(
            gene,
            reference_haplotype_name,
            haplotypes,
            frozenset([rs_id_info4]),
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                haplotypes,
                frozenset([rs_id_info1, rs_id_info2]),
                drugs,
                rs_id_to_ref_seq_diff_annotation,
            )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                haplotypes,
                frozenset([rs_id_info1, rs_id_info3]),
                drugs,
                rs_id_to_ref_seq_diff_annotation,
            )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                haplotypes,
                frozenset([rs_id_info1, rs_id_info2, rs_id_info3, rs_id_info4]),
                drugs,
                rs_id_to_ref_seq_diff_annotation,
            )

    def test_gene_info_with_rs_id_in_haplotype_without_info(self) -> None:
        """Error when haplotype uses rs id for which there is no RsIdInfo object"""
        gene = "FAKE"
        reference_haplotype_name = "*1"
        rs_id_infos: FrozenSet[RsIdInfo] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_diff_annotation: Dict[str, Annotation] = dict()

        empty_haplotypes: FrozenSet[Haplotype] = frozenset()
        non_empty_haplotypes = frozenset([Haplotype("*2", "No Function", frozenset([Variant("rs238423", "A")]))])

        GeneInfo(gene, reference_haplotype_name, empty_haplotypes, rs_id_infos, drugs, rs_id_to_ref_seq_diff_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                non_empty_haplotypes,
                rs_id_infos,
                drugs,
                rs_id_to_ref_seq_diff_annotation,
            )

    def test_gene_info_with_variant_that_is_ref_v38(self) -> None:
        """Error when gene info has variant where variant allele is the ref allele"""
        gene = "FAKE"
        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        drugs: FrozenSet[DrugInfo] = frozenset()

        rs_id = "rs294924"

        rs_id_to_ref_seq_diff_annotation = {rs_id: Annotation("399483A>C", "399483C>A")}
        haplotypes = frozenset([Haplotype("*3", "No Function", frozenset([Variant("rs294924", "G")]))])

        rs_id_infos1 = frozenset(
            [
                RsIdInfo(
                    rs_id,
                    ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "A"),
                    ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "C"),
                )
            ]
        )
        rs_id_infos2 = frozenset(
            [
                RsIdInfo(
                    rs_id,
                    ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "A"),
                    ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "G"),
                )
            ]
        )

        GeneInfo(gene, reference_haplotype_name, haplotypes, rs_id_infos1, drugs, rs_id_to_ref_seq_diff_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, reference_haplotype_name, haplotypes, rs_id_infos2, drugs, rs_id_to_ref_seq_diff_annotation)

    def test_gene_info_with_ref_seq_difference_without_annotation(self) -> None:
        """Error when a ref seq difference does not have an annotation"""
        gene = "FAKE"
        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_diff_annotation: Dict[str, Annotation] = dict()
        rs_id_info = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "A"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "C"),
        )

        empty_rs_id_infos: FrozenSet[RsIdInfo] = frozenset()
        non_empty_rs_id_infos = frozenset([rs_id_info])

        GeneInfo(gene, reference_haplotype_name, haplotypes, empty_rs_id_infos, drugs, rs_id_to_ref_seq_diff_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                haplotypes,
                non_empty_rs_id_infos,
                drugs,
                rs_id_to_ref_seq_diff_annotation,
            )

    def test_gene_info_with_extra_ref_seq_difference_annotation(self) -> None:
        """Error when a ref seq difference annotation does not match known a known ref seq difference"""
        gene = "FAKE"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        rs_id_infos: FrozenSet[RsIdInfo] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()

        empty_rs_id_to_ref_seq_diff_annotation: Dict[str, Annotation] = dict()
        non_empty_rs_id_to_ref_seq_diff_annotation = {"rs2493023": Annotation("3445T>A", "3445A>T")}

        GeneInfo(gene, reference_haplotype_name, haplotypes, rs_id_infos, drugs, empty_rs_id_to_ref_seq_diff_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                haplotypes,
                rs_id_infos,
                drugs,
                non_empty_rs_id_to_ref_seq_diff_annotation,
            )

    def test_haplotype_with_repeat_rs_ids(self) -> None:
        """Error when haplotype has multiple variants with the same rs id"""
        name = "*5"
        function = "No Function"

        non_clashing_variants = frozenset([Variant("rs88293", "A"), Variant("rs39492", "T")])
        clashing_variants = frozenset([Variant("rs88293", "A"), Variant("rs88293", "T")])

        Haplotype(name, function, non_clashing_variants)
        with self.assertRaises(ValueError):
            Haplotype(name, function, clashing_variants)

    def test_empty_haplotype(self) -> None:
        """Error when haplotype has no variants"""
        name = "*5"
        function = "No Function"

        variants = frozenset([Variant("rs88293", "A"), Variant("rs39492", "T")])
        empty_variants: FrozenSet[Variant] = frozenset()

        Haplotype(name, function, variants)
        with self.assertRaises(ValueError):
            Haplotype(name, function, empty_variants)

    def test_panel_with_repeat_gene_names(self) -> None:
        """Error when panel has info for gene multiple times"""
        name = "FakePanel"
        version = "1.0"

        gene1 = "FAKE"
        gene2 = "OTHER"

        reference_haplotype_name = "*1"
        other_reference_haplotype_name = "*1_something else"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        rs_id_infos: FrozenSet[RsIdInfo] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_diff_annotation: Dict[str, Annotation] = dict()

        gene_info1 = GeneInfo(
            gene1,
            reference_haplotype_name,
            haplotypes,
            rs_id_infos,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        gene_info2 = GeneInfo(
            gene2,
            reference_haplotype_name,
            haplotypes,
            rs_id_infos,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        gene_info3 = GeneInfo(
            gene1,
            other_reference_haplotype_name,
            haplotypes,
            rs_id_infos,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        gene_infos_without_repeat = frozenset([gene_info1, gene_info2])
        gene_infos_with_repeat = frozenset([gene_info1, gene_info3])

        Panel(name, version, gene_infos_without_repeat)
        with self.assertRaises(ValueError):
            Panel(name, version, gene_infos_with_repeat)

    def test_panel_with_overlapping_rs_id_infos_for_different_genes(self) -> None:
        """Error when panel has overlapping rs id infos for different genes, but not when they are exactly the same"""
        name = "FakePanel"
        version = "1.0"

        gene1 = "FAKE"
        gene2 = "OTHER"

        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_diff_annotation: Dict[str, Annotation] = dict()

        rs_id_info1 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "AT"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "AT"),
        )
        rs_id_info2 = RsIdInfo(
            "rs3949923",
            ReferenceSite(GeneCoordinate(chromosome_v37, 293993), "C"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 1388323), "C"),
        )
        rs_id_info3 = RsIdInfo(
            "rs12993",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499592), "GG"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399482), "GG"),
        )

        rs_id_infos1 = frozenset([rs_id_info1])
        rs_id_infos2 = frozenset([rs_id_info1, rs_id_info2])
        rs_id_infos3 = frozenset([rs_id_info3])
        gene_info1 = GeneInfo(
            gene1,
            reference_haplotype_name,
            haplotypes,
            rs_id_infos1,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        gene_info2 = GeneInfo(
            gene2,
            reference_haplotype_name,
            haplotypes,
            rs_id_infos2,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        gene_info3 = GeneInfo(
            gene2,
            reference_haplotype_name,
            haplotypes,
            rs_id_infos3,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )

        Panel(name, version, frozenset([gene_info1, gene_info2]))
        with self.assertRaises(ValueError):
            Panel(name, version, frozenset([gene_info1, gene_info3]))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
