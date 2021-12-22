import unittest
from typing import Dict, FrozenSet

from base.gene_coordinate import GeneCoordinate
from base.reference_site import ReferenceSite
from panel.annotation import Annotation
from panel.drug_info import DrugInfo
from panel.gene_info import GeneInfo
from panel.haplotype import Haplotype
from panel.panel import Panel
from panel.rs_id_info import RsIdInfo
from panel.variant import Variant


class TestConfig(unittest.TestCase):
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

        GeneInfo(
            gene,
            reference_haplotype_name,
            None,
            haplotypes1,
            rs_id_infos,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                None,
                haplotypes2,
                rs_id_infos,
                drugs,
                rs_id_to_ref_seq_diff_annotation,
            )

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

        GeneInfo(
            gene,
            reference_haplotype_name,
            None,
            haplotypes1,
            rs_id_infos,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                None,
                haplotypes2,
                rs_id_infos,
                drugs,
                rs_id_to_ref_seq_diff_annotation,
            )

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

        GeneInfo(gene, reference_haplotype_name, None, haplotypes, rs_id_infos, drugs, rs_id_to_ref_seq_diff_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                None,
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

        GeneInfo(
            gene,
            reference_haplotype_name,
            None,
            haplotypes,
            single_rs_id_info,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                None,
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
            None,
            haplotypes,
            frozenset([rs_id_info1]),
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        GeneInfo(
            gene,
            reference_haplotype_name,
            None,
            haplotypes,
            frozenset([rs_id_info2]),
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        GeneInfo(
            gene,
            reference_haplotype_name,
            None,
            haplotypes,
            frozenset([rs_id_info3]),
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        GeneInfo(
            gene,
            reference_haplotype_name,
            None,
            haplotypes,
            frozenset([rs_id_info4]),
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                None,
                haplotypes,
                frozenset([rs_id_info1, rs_id_info2]),
                drugs,
                rs_id_to_ref_seq_diff_annotation,
            )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                None,
                haplotypes,
                frozenset([rs_id_info1, rs_id_info3]),
                drugs,
                rs_id_to_ref_seq_diff_annotation,
            )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                None,
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

        GeneInfo(
            gene,
            reference_haplotype_name,
            None,
            empty_haplotypes,
            rs_id_infos,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                None,
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

        GeneInfo(
            gene,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos1,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                None,
                haplotypes,
                rs_id_infos2,
                drugs,
                rs_id_to_ref_seq_diff_annotation,
            )

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

        GeneInfo(
            gene,
            reference_haplotype_name,
            "ENST000000384893",
            haplotypes,
            empty_rs_id_infos,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                "ENST000000384893",
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

        GeneInfo(
            gene,
            reference_haplotype_name,
            "ENST000000384893",
            haplotypes,
            rs_id_infos,
            drugs,
            empty_rs_id_to_ref_seq_diff_annotation,
        )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                "ENST000000384893",
                haplotypes,
                rs_id_infos,
                drugs,
                non_empty_rs_id_to_ref_seq_diff_annotation,
            )

    def test_gene_info_with_repeated_rs_ids(self) -> None:
        """Error when a gene info contains an rs id multiple times"""
        gene = "FAKE"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_diff_annotation: Dict[str, Annotation] = dict()

        chromosome_v37 = "X"
        chromosome_v38 = "chrX"
        rs_id_info1 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "AT"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "AT"),
        )
        rs_id_info2 = RsIdInfo(
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499594), "AT"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "AT"),
        )

        rs_id_infos1 = frozenset([rs_id_info1])
        rs_id_infos2 = frozenset([rs_id_info1, rs_id_info2])

        GeneInfo(
            gene,
            reference_haplotype_name,
            "ENST000000384893",
            haplotypes,
            rs_id_infos1,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        with self.assertRaises(ValueError):
            GeneInfo(
                gene,
                reference_haplotype_name,
                "ENST000000384893",
                haplotypes,
                rs_id_infos2,
                drugs,
                rs_id_to_ref_seq_diff_annotation,
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
            None,
            haplotypes,
            rs_id_infos,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        gene_info2 = GeneInfo(
            gene2,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        gene_info3 = GeneInfo(
            gene1,
            other_reference_haplotype_name,
            None,
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
        """Error when panel has overlapping rs id infos for different genes"""
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
        rs_id_infos2 = frozenset([rs_id_info2])
        rs_id_infos3 = frozenset([rs_id_info3])
        gene_info1 = GeneInfo(
            gene1,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos1,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        gene_info2 = GeneInfo(
            gene2,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos2,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        gene_info3 = GeneInfo(
            gene2,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos3,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )

        Panel(name, version, frozenset([gene_info1, gene_info2]))
        with self.assertRaises(ValueError):
            Panel(name, version, frozenset([gene_info1, gene_info3]))

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
            "rs294924",
            ReferenceSite(GeneCoordinate(chromosome_v37, 499593), "AT"),
            ReferenceSite(GeneCoordinate(chromosome_v38, 399483), "AT"),
        )

        rs_id_infos1 = frozenset([rs_id_info1])
        rs_id_infos2 = frozenset([rs_id_info2])
        rs_id_infos3 = frozenset([rs_id_info3])
        gene_info1 = GeneInfo(
            gene1,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos1,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        gene_info2 = GeneInfo(
            gene2,
            reference_haplotype_name,
            None,
            haplotypes,
            rs_id_infos2,
            drugs,
            rs_id_to_ref_seq_diff_annotation,
        )
        gene_info3 = GeneInfo(
            gene2,
            reference_haplotype_name,
            None,
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
