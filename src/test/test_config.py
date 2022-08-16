import unittest
from typing import FrozenSet

from util.gene_coordinate import GeneCoordinate
from util.reference_site import ReferenceSite
from panel.annotation import Annotation
from panel.drug_summary import DrugSummary
from panel.gene_panel import GenePanel
from panel.haplotype import Haplotype
from panel.panel import Panel
from panel.rs_id_info import RsIdInfo
from panel.variant import Variant


class TestConfig(unittest.TestCase):
    def test_rs_id_info_with_extra_ref_seq_difference_annotation(self) -> None:
        """Error when an rs id has a ref sequence difference annotation set when the alleles are identical"""
        rs_id = "rs2493023"
        reference_site_v37 = ReferenceSite(GeneCoordinate("1", 3758), "C")
        reference_site_v38 = ReferenceSite(GeneCoordinate("1", 3993), "C")
        ref_seq_difference_annotation = Annotation("3758C>A", "3993C>A")
        RsIdInfo(
            rs_id,
            reference_site_v37,
            reference_site_v38,
            None,
        )
        with self.assertRaises(ValueError):
            RsIdInfo(
                rs_id,
                reference_site_v37,
                reference_site_v38,
                ref_seq_difference_annotation,
            )

    def test_rs_id_info_with_missing_ref_seq_difference_annotation(self) -> None:
        """Error when an rs id does not have a ref sequence difference annotation set when the alleles differ"""
        rs_id = "rs2493023"
        reference_site_v37 = ReferenceSite(GeneCoordinate("1", 3758), "C")
        reference_site_v38 = ReferenceSite(GeneCoordinate("1", 3993), "A")
        ref_seq_difference_annotation = Annotation("3758C>A", "3993A>C")
        RsIdInfo(
            rs_id,
            reference_site_v37,
            reference_site_v38,
            ref_seq_difference_annotation,
        )
        with self.assertRaises(ValueError):
            RsIdInfo(
                rs_id,
                reference_site_v37,
                reference_site_v38,
                None,
            )

    def test_gene_panel_with_overlapping_haplotype_names(self) -> None:
        """Error when haplotype name used multiple times for gene"""
        gene = "FAKE"
        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        drugs: FrozenSet[DrugSummary] = frozenset()

        variant1 = Variant("rs94982", "A")
        variant2 = Variant("rs394934", "T")

        variant1_rs_id_info = RsIdInfo(
            variant1.rs_id,
            ReferenceSite(GeneCoordinate(chromosome_v37, 4994545), "C"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 2993823), "C"),
            None,
        )
        variant2_rs_id_info = RsIdInfo(
            variant2.rs_id,
            ReferenceSite(GeneCoordinate(chromosome_v37, 3993842), "G"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 2949923), "G"),
            None,
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

        GenePanel(
            gene,
            reference_haplotype_name,
            None,
            haplotypes1,
            rs_id_infos,
            drugs,
        )
        with self.assertRaises(ValueError):
            GenePanel(
                gene,
                reference_haplotype_name,
                None,
                haplotypes2,
                rs_id_infos,
                drugs,
            )

    def test_gene_panel_with_overlapping_haplotype_variants(self) -> None:
        """Error when different haplotypes have the exact same variant combination"""
        gene = "FAKE"
        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        drugs: FrozenSet[DrugSummary] = frozenset()

        variant1 = Variant("rs94982", "A")
        variant2 = Variant("rs394934", "T")
        variant3 = Variant("rs495825", "C")

        variant1_rs_id_info = RsIdInfo(
            variant1.rs_id,
            ReferenceSite(GeneCoordinate(chromosome_v37, 4994545), "C"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 2993823), "C"),
            None,
        )
        variant2_rs_id_info = RsIdInfo(
            variant2.rs_id,
            ReferenceSite(GeneCoordinate(chromosome_v37, 3993842), "G"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 2949923), "G"),
            None,
        )
        variant3_rs_id_info = RsIdInfo(
            variant3.rs_id,
            ReferenceSite(GeneCoordinate(chromosome_v37, 293923), "A"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 138812), "A"),
            None,
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

        GenePanel(
            gene,
            reference_haplotype_name,
            None,
            haplotypes1,
            rs_id_infos,
            drugs,
        )
        with self.assertRaises(ValueError):
            GenePanel(
                gene,
                reference_haplotype_name,
                None,
                haplotypes2,
                rs_id_infos,
                drugs,
            )

    def test_gene_panel_with_overlapping_drug_names(self) -> None:
        """Error when drug name used multiple times for gene"""
        gene = "FAKE"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        rs_id_infos: FrozenSet[RsIdInfo] = frozenset()

        drug_summary1 = DrugSummary("Paracetamol", "google.com")
        drug_summary2 = DrugSummary("Paracetamol", "wikipedia.com")
        drugs = frozenset([drug_summary1])
        overlapping_drugs = frozenset([drug_summary1, drug_summary2])

        GenePanel(gene, reference_haplotype_name, None, haplotypes, rs_id_infos, drugs)
        with self.assertRaises(ValueError):
            GenePanel(
                gene,
                reference_haplotype_name,
                None,
                haplotypes,
                rs_id_infos,
                overlapping_drugs,
            )

    def test_gene_panel_with_overlapping_rs_id_infos(self) -> None:
        """Error when gene info has rs id infos for which the relevant coordinates overlap"""
        gene = "FAKE"
        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugSummary] = frozenset()

        rs_id_info1 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "A"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "A"),
            None,
        )
        rs_id_info2 = RsIdInfo(
            "rs294927",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499592), "AA"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399482), "AA"),
            None,
        )
        single_rs_id_info = frozenset([rs_id_info1])
        overlapping_rs_id_infos = frozenset([rs_id_info1, rs_id_info2])

        GenePanel(
            gene,
            reference_haplotype_name,
            None,
            haplotypes,
            single_rs_id_info,
            drugs,
        )
        with self.assertRaises(ValueError):
            GenePanel(
                gene,
                reference_haplotype_name,
                None,
                haplotypes,
                overlapping_rs_id_infos,
                drugs,
            )

    def test_gene_panel_with_rs_id_infos_for_different_chromosome(self) -> None:
        """Error when gene info has rs id infos for which the relevant coordinates have a different chromosome"""
        gene = "FAKE"
        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugSummary] = frozenset()

        other_chromosome_v37 = "1"
        other_chromosome_v38 = "1"

        rs_id_info1 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "A"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "A"),
            None,
        )
        rs_id_info2 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(other_chromosome_v37, 499593), "A"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "A"),
            None,
        )
        rs_id_info3 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "A"),
            ReferenceSite(GeneCoordinate(other_chromosome_v38, 399483), "A"),
            None,
        )
        rs_id_info4 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(other_chromosome_v37, 499593), "A"),
            ReferenceSite(GeneCoordinate(other_chromosome_v38, 399483), "A"),
            None,
        )

        GenePanel(
            gene,
            reference_haplotype_name,
            None,
            haplotypes,
            frozenset([rs_id_info1]),
            drugs,
        )
        GenePanel(
            gene,
            reference_haplotype_name,
            None,
            haplotypes,
            frozenset([rs_id_info2]),
            drugs,
        )
        GenePanel(
            gene,
            reference_haplotype_name,
            None,
            haplotypes,
            frozenset([rs_id_info3]),
            drugs,
        )
        GenePanel(
            gene,
            reference_haplotype_name,
            None,
            haplotypes,
            frozenset([rs_id_info4]),
            drugs,
        )
        with self.assertRaises(ValueError):
            GenePanel(
                gene,
                reference_haplotype_name,
                None,
                haplotypes,
                frozenset([rs_id_info1, rs_id_info2]),
                drugs,
            )
        with self.assertRaises(ValueError):
            GenePanel(
                gene,
                reference_haplotype_name,
                None,
                haplotypes,
                frozenset([rs_id_info1, rs_id_info3]),
                drugs,
            )
        with self.assertRaises(ValueError):
            GenePanel(
                gene,
                reference_haplotype_name,
                None,
                haplotypes,
                frozenset([rs_id_info1, rs_id_info2, rs_id_info3, rs_id_info4]),
                drugs,
            )

    def test_gene_panel_with_rs_id_in_haplotype_without_info(self) -> None:
        """Error when haplotype uses rs id for which there is no RsIdInfo object"""
        gene = "FAKE"
        reference_haplotype_name = "*1"
        rs_id_infos: FrozenSet[RsIdInfo] = frozenset()
        drugs: FrozenSet[DrugSummary] = frozenset()

        empty_haplotypes: FrozenSet[Haplotype] = frozenset()
        non_empty_haplotypes = frozenset([Haplotype("*2", "No Function", frozenset([Variant("rs238423", "A")]))])

        GenePanel(
            gene,
            reference_haplotype_name,
            None,
            empty_haplotypes,
            rs_id_infos,
            drugs,
        )
        with self.assertRaises(ValueError):
            GenePanel(
                gene,
                reference_haplotype_name,
                None,
                non_empty_haplotypes,
                rs_id_infos,
                drugs,
            )

    def test_gene_panel_with_variant_that_is_ref_v38(self) -> None:
        """Error when gene info has variant where variant allele is the ref allele"""
        gene = "FAKE"
        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        drugs: FrozenSet[DrugSummary] = frozenset()

        rs_id = "rs294924"

        haplotypes = frozenset([Haplotype("*3", "No Function", frozenset([Variant("rs294924", "G")]))])

        rs_id_infos1 = frozenset(
            [
                RsIdInfo(
                    rs_id,
                    ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "A"),
                    ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "C"),
                    Annotation("399483A>C", "399483C>A")
                )
            ]
        )
        rs_id_infos2 = frozenset(
            [
                RsIdInfo(
                    rs_id,
                    ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "A"),
                    ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "G"),
                    Annotation("399483A>C", "399483C>A")
                )
            ]
        )

        GenePanel(
            gene,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos1,
            drugs,
        )
        with self.assertRaises(ValueError):
            GenePanel(
                gene,
                reference_haplotype_name,
                None,
                haplotypes,
                rs_id_infos2,
                drugs,
            )

    def test_gene_panel_with_repeated_rs_ids(self) -> None:
        """Error when a gene info contains an rs id multiple times"""
        gene = "FAKE"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugSummary] = frozenset()

        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        rs_id_info1 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "AT"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "AT"),
            None,
        )
        rs_id_info2 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499594), "AT"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "AT"),
            None,
        )

        rs_id_infos1 = frozenset([rs_id_info1])
        rs_id_infos2 = frozenset([rs_id_info1, rs_id_info2])

        GenePanel(
            gene,
            reference_haplotype_name,
            "ENST000000384893",
            haplotypes,
            rs_id_infos1,
            drugs,
        )
        with self.assertRaises(ValueError):
            GenePanel(
                gene,
                reference_haplotype_name,
                "ENST000000384893",
                haplotypes,
                rs_id_infos2,
                drugs,
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
        drugs: FrozenSet[DrugSummary] = frozenset()

        gene_panel1 = GenePanel(
            gene1,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos,
            drugs,
        )
        gene_panel2 = GenePanel(
            gene2,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos,
            drugs,
        )
        gene_panel3 = GenePanel(
            gene1,
            other_reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos,
            drugs,
        )
        gene_panels_without_repeat = frozenset([gene_panel1, gene_panel2])
        gene_panels_with_repeat = frozenset([gene_panel1, gene_panel3])

        Panel(name, version, gene_panels_without_repeat)
        with self.assertRaises(ValueError):
            Panel(name, version, gene_panels_with_repeat)

    def test_panel_with_overlapping_rs_id_infos_for_different_genes(self) -> None:
        """Error when panel has overlapping rs id infos for different genes"""
        name = "FakePanel"
        version = "1.0"

        gene1 = "FAKE"
        gene2 = "OTHER"

        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugSummary] = frozenset()

        rs_id_info1 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "AT"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "AT"),
            None,
        )
        rs_id_info2 = RsIdInfo(
            "rs3949923",
            ReferenceSite(GeneCoordinate(chromosome_v37, 293993), "C"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 1388323), "C"),
            None,
        )
        rs_id_info3 = RsIdInfo(
            "rs12993",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499592), "GG"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399482), "GG"),
            None,
        )

        rs_id_infos1 = frozenset([rs_id_info1])
        rs_id_infos2 = frozenset([rs_id_info2])
        rs_id_infos3 = frozenset([rs_id_info3])
        gene_panel1 = GenePanel(
            gene1,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos1,
            drugs,
        )
        gene_panel2 = GenePanel(
            gene2,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos2,
            drugs,
        )
        gene_panel3 = GenePanel(
            gene2,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos3,
            drugs,
        )

        Panel(name, version, frozenset([gene_panel1, gene_panel2]))
        with self.assertRaises(ValueError):
            Panel(name, version, frozenset([gene_panel1, gene_panel3]))

    def test_panel_with_repeated_rs_ids(self) -> None:
        """Error when panel has repeated rs ids for different genes"""
        name = "FakePanel"
        version = "1.0"

        gene1 = "FAKE"
        gene2 = "OTHER"

        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugSummary] = frozenset()

        rs_id_info1 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "AT"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "AT"),
            None,
        )
        rs_id_info2 = RsIdInfo(
            "rs3949923",
            ReferenceSite(GeneCoordinate(chromosome_v37, 293993), "C"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 1388323), "C"),
            None,
        )
        rs_id_info3 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "AT"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "AT"),
            None,
        )

        rs_id_infos1 = frozenset([rs_id_info1])
        rs_id_infos2 = frozenset([rs_id_info2])
        rs_id_infos3 = frozenset([rs_id_info3])
        gene_panel1 = GenePanel(
            gene1,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos1,
            drugs,
        )
        gene_panel2 = GenePanel(
            gene2,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos2,
            drugs,
        )
        gene_panel3 = GenePanel(
            gene2,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos3,
            drugs,
        )

        Panel(name, version, frozenset([gene_panel1, gene_panel2]))
        with self.assertRaises(ValueError):
            Panel(name, version, frozenset([gene_panel1, gene_panel3]))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
