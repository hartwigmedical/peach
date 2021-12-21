import unittest
from typing import Dict, Set

from base.filter import FullCallFilter, SimpleCallFilter
from base.gene_coordinate import GeneCoordinate
from base.reference_assembly import ReferenceAssembly
from base.reference_site import ReferenceSite
from call_data import SimpleCallData, SimpleCall, HaplotypeCall, FullCall, FullCallData
from config.annotation import Annotation
from config.drug_info import DrugInfo
from config.gene_info import GeneInfo
from config.haplotype import Haplotype
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant
from analysis.pgx_analysis import PgxAnalyser, PgxAnalysis


class TestPgxAnalysis(unittest.TestCase):
    @classmethod
    def __get_wide_example_panel(cls) -> Panel:
        dpyd_two_a_variant = Variant("rs3918290", "T")
        dpyd_two_b_variant = Variant("rs1801159", "C")
        dpyd_three_variant = Variant("rs72549303", "TG")
        fake_variant = Variant("rs1212125", "C")
        fake2_variant = Variant("rs1212127", "C")

        dpyd_haplotypes = frozenset({
            Haplotype("*2A", "No Function", frozenset({dpyd_two_a_variant})),
            Haplotype("*2B", "No Function", frozenset({dpyd_two_a_variant, dpyd_two_b_variant})),
            Haplotype("*3", "Normal Function", frozenset({dpyd_three_variant})),
        })
        dpyd_rs_id_infos = frozenset({
            RsIdInfo("rs3918290", ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C")),
            RsIdInfo("rs72549309", ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA")),
            RsIdInfo("rs1801159", ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T")),
            RsIdInfo("rs72549303", ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC")),
        })
        dpyd_drugs = frozenset({
            DrugInfo("5-Fluorouracil", "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"),
            DrugInfo("Capecitabine", "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"),
        })
        dpyd_rs_id_to_ref_seq_difference_annotations = {
            "rs72549303": Annotation("6744CA>GA", "6744GA>CA"),
        }

        fake_haplotypes = frozenset({
            Haplotype("*4A", "Reduced Function", frozenset({fake_variant})),
        })
        fake_rs_id_infos = frozenset({
            RsIdInfo("rs1212125", ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T")),
        })
        fake_drugs = frozenset({
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
        })
        fake_rs_id_to_ref_seq_difference_annotations: Dict[str, Annotation] = {}

        fake2_haplotypes = frozenset({
            Haplotype("*4A", "Reduced Function", frozenset({fake2_variant})),
        })
        fake2_rs_id_infos = frozenset({
            RsIdInfo("rs1212127", ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T")),
        })
        fake2_drugs = frozenset({
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
        })
        fake2_rs_id_to_ref_seq_difference_annotations = {"rs1212127": Annotation("1324C>T", "1324T>C")}

        dpyd_gene_info = GeneInfo(
            "DPYD", "*1", "ENST00000370192", dpyd_haplotypes, dpyd_rs_id_infos,
            dpyd_drugs, dpyd_rs_id_to_ref_seq_difference_annotations,
        )
        fake_gene_info = GeneInfo(
            "FAKE", "*1", None, fake_haplotypes, fake_rs_id_infos,
            fake_drugs, fake_rs_id_to_ref_seq_difference_annotations,
        )
        fake2_gene_info = GeneInfo(
            "FAKE2", "*1", None, fake2_haplotypes, fake2_rs_id_infos,
            fake2_drugs, fake2_rs_id_to_ref_seq_difference_annotations,
        )
        gene_infos = frozenset({dpyd_gene_info, fake_gene_info, fake2_gene_info})

        name = "WideTestPanel"
        version = "1.0"
        return Panel(name, version, gene_infos)

    @classmethod
    def __get_narrow_example_panel(cls, included_haplotypes: Set[str]) -> Panel:
        dpyd_two_a_variant = Variant("rs3918290", "T")
        dpyd_five_variant = Variant("rs1801159", "C")
        dpyd_three_variant = Variant("rs72549303", "TG")
        dpyd_seven_variant = Variant("rs2938101", "AGT")

        possible_dpyd_haplotypes = [
            Haplotype("*2A", "No Function", frozenset({dpyd_two_a_variant})),
            Haplotype("*2B", "No Function", frozenset({dpyd_two_a_variant, dpyd_five_variant})),
            Haplotype("*3", "Normal Function", frozenset({dpyd_three_variant})),
            Haplotype("*5", "Normal Function", frozenset({dpyd_five_variant})),
            Haplotype("*7", "Normal Function", frozenset({dpyd_seven_variant})),
            Haplotype("*9", "Normal Function", frozenset({dpyd_two_a_variant, dpyd_three_variant})),
            Haplotype("*10", "Normal Function", frozenset({dpyd_three_variant, dpyd_five_variant})),
        ]
        possible_rs_id_infos = [
            RsIdInfo("rs3918290", ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C")),
            RsIdInfo("rs72549309", ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA")),
            RsIdInfo("rs1801159", ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T")),
            RsIdInfo("rs72549303", ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC")),
            RsIdInfo("rs2938101", ReferenceSite(GeneCoordinate("1", 97912838), "A"), ReferenceSite(GeneCoordinate("chr1", 97453984), "A")),
        ]

        unknown_haplotypes = included_haplotypes.difference({haplotype.name for haplotype in possible_dpyd_haplotypes})
        if unknown_haplotypes:
            raise ValueError(f"Not all dpyd haplotype names recognized: unknown={unknown_haplotypes}")

        # use a set to cause the order of haplotypes to be random, to make sure this order will not matter.
        included_dpyd_haplotypes = frozenset(
            {haplotype for haplotype in possible_dpyd_haplotypes if haplotype.name in included_haplotypes}
        )
        included_dpyd_rs_ids = {"rs72549303"}.union(
            {variant.rs_id for haplotype in included_dpyd_haplotypes for variant in haplotype.variants}
        )
        included_dpyd_rs_id_infos = frozenset(
            {rs_id_info for rs_id_info in possible_rs_id_infos if rs_id_info.rs_id in included_dpyd_rs_ids}
        )
        dpyd_drugs = frozenset({
            DrugInfo("5-Fluorouracil", "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"),
            DrugInfo("Capecitabine", "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"),
        })
        dpyd_rs_id_to_difference_annotations = {
            "rs72549303": Annotation("6744CA>GA", "6744GA>CA"),
        }

        gene_infos = frozenset({
            GeneInfo(
                "DPYD", "*1", "ENST00000370192", included_dpyd_haplotypes, included_dpyd_rs_id_infos,
                dpyd_drugs, dpyd_rs_id_to_difference_annotations,
            ),
        })
        name = "NarrowTestPanel"
        version = "1.0"
        return Panel(name, version, gene_infos)

    def test_empty_v37(self) -> None:
        """No variants wrt v37"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset(), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*3", 2)}, "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": {HaplotypeCall("*4A", 2)}
        }

        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TG"), "DPYD", ("rs72549303",),
                "REF_CALL", FullCallFilter.NO_CALL, "6744GA>CA", FullCallFilter.INFERRED_PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("C", "C"), "FAKE2", ("rs1212127",),
                "REF_CALL", FullCallFilter.NO_CALL, "1324T>C", FullCallFilter.INFERRED_PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "T"), "FAKE", ("rs1212125",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_empty_v38(self) -> None:
        """No variants wrt v38"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset(), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*1", 2)}, "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": {HaplotypeCall("*1", 2)}
        }

        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.INFERRED_PASS, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("T", "T"), "FAKE2", ("rs1212127",),
                "1324C>T", FullCallFilter.INFERRED_PASS, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "T"), "FAKE", ("rs1212125",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_hom_ref_from_v37(self) -> None:
        """All haplotypes are  *1_HOM with v37 input data"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("T", "T"),
                "FAKE2", ("rs1212127",), "1324C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TC", "TC"),
                "DPYD", (".",), "6744CA>GA", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("C", "C"),
                "DPYD", ("rs3918290",), "REF_CALL", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*1", 2)}, "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": {HaplotypeCall("*1", 2)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("T", "T"), "FAKE2", ("rs1212127",),
                "1324C>T", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "T"), "FAKE", ("rs1212125",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_heterozygous_v37(self) -> None:
        """All haplotypes are heterozygous, both with single and multiple non-ref haplotypes. Uses v37 input"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("C", "T"),
                "FAKE2", ("rs1212127",), "1324C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("C", "T"),
                "DPYD", ("rs3918290",), "35G>A", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*2B", 1), HaplotypeCall("*3", 1)},
            "FAKE": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)},
            "FAKE2": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)},
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "T"), "DPYD", ("rs3918290",),
                "35G>A", FullCallFilter.PASS, "35G>A", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "6744GA>CA", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "C"), "DPYD", ("rs1801159",),
                "674A>G", FullCallFilter.PASS, "674A>G", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("C", "T"), "FAKE2", ("rs1212127",),
                "1324C>T", FullCallFilter.PASS, "1324T>C", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "C"), "FAKE", ("rs1212125",),
                "1005T>C", FullCallFilter.PASS, "1005T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_heterozygous_v38(self) -> None:
        """All haplotypes are heterozygous, both with single and multiple non-ref haplotypes. Uses v38 input"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "T"),
                "FAKE2", ("rs1212127",), "1324T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr5", 97450060), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744GA>CA", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("C", "T"),
                "DPYD", ("rs3918290",), "35G>A", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*2B", 1), HaplotypeCall("*3", 1)},
            "FAKE": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)},
            "FAKE2": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)},
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "T"), "DPYD", ("rs3918290",),
                "35G>A", FullCallFilter.PASS, "35G>A", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "6744GA>CA", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "C"), "DPYD", ("rs1801159",),
                "674A>G", FullCallFilter.PASS, "674A>G", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("C", "T"), "FAKE2", ("rs1212127",),
                "1324C>T", FullCallFilter.PASS, "1324T>C", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "C"), "FAKE", ("rs1212125",),
                "1005T>C", FullCallFilter.PASS, "1005T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_more_than_two_haplotypes_present(self) -> None:
        """More than two haplotypes are present, specifically three. Uses v37 input"""
        panel = self.__get_narrow_example_panel({"*2A", "*3", "*7"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("C", "T"),
                "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97912838), "A"), ("AGT", "AGT"),
                "DPYD", ("rs2938101",), "293A>AGT", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*2A", 1), HaplotypeCall("*3", 2), HaplotypeCall("*7", 2)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97912838), "A"), ReferenceSite(GeneCoordinate("chr1", 97453984), "A"),
                ("AGT", "AGT"), "DPYD", ("rs2938101",),
                "293A>AGT", FullCallFilter.PASS, "293A>AGT", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "T"), "DPYD", ("rs3918290",),
                "9213C>T", FullCallFilter.PASS, "9213C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TG"), "DPYD", ("rs72549303",),
                "REF_CALL", FullCallFilter.NO_CALL, "6744GA>CA", FullCallFilter.INFERRED_PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ambiguous_haplotype_with_clear_winner_homozygous(self) -> None:
        """Ambiguous homozygous haplotype, where the simpler possibility should be preferred. Uses v38 input"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "REF_CALL", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {"DPYD": {HaplotypeCall("*2B", 2)}}
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("T", "T"), "DPYD", ("rs3918290",),
                "9213C>T", FullCallFilter.PASS, "9213C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("C", "C"), "DPYD", ("rs1801159",),
                "293T>C", FullCallFilter.PASS, "293T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ambiguous_haplotype_with_clear_winner_heterozygous(self) -> None:
        """Ambiguous heterozygous haplotype, where the simpler possibility should be preferred. Uses v37 input."""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("C", "T"),
                "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "293T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {"DPYD": {HaplotypeCall("*2B", 1), HaplotypeCall("*1", 1)}}
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "T"), "DPYD", ("rs3918290",),
                "9213C>T", FullCallFilter.PASS, "9213C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "C"), "DPYD", ("rs1801159",),
                "293T>C", FullCallFilter.PASS, "293T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ambiguous_haplotype_with_clear_winner_mix(self) -> None:
        """
        Ambiguous mix of homozygous and heterozygous haplotype, where the simplest possibility should be preferred.
        Uses v38 input
        """
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "293T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "REF_CALL", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {"DPYD": {HaplotypeCall("*2B", 1), HaplotypeCall("*2A", 1)}}
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("T", "T"), "DPYD", ("rs3918290",),
                "9213C>T", FullCallFilter.PASS, "9213C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "C"), "DPYD", ("rs1801159",),
                "293T>C", FullCallFilter.PASS, "293T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_complicated_ambiguous_haplotype_with_a_clear_winner(self) -> None:
        """Ambiguous set of haplotypes, with no haplotype that is simplest. Uses v37 input"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10", "*7"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "293T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97912838), "A"), ("AGT", "AGT"),
                "DPYD", ("rs2938101",), "301A>AGT", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*2B", 1), HaplotypeCall("*9", 1), HaplotypeCall("*7", 2)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97912838), "A"), ReferenceSite(GeneCoordinate("chr1", 97453984), "A"),
                ("AGT", "AGT"), "DPYD", ("rs2938101",),
                "301A>AGT", FullCallFilter.PASS, "301A>AGT", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("T", "T"), "DPYD", ("rs3918290",),
                "9213C>T", FullCallFilter.PASS, "9213C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "6744GA>CA", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "C"), "DPYD", ("rs1801159",),
                "293T>C", FullCallFilter.PASS, "293T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ambiguous_homozygous_haplotype_with_a_less_clear_winner(self) -> None:
        """Ambiguous set of homozygous haplotypes, with winning haplotype that might not be ideal. Uses v38 input"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TG"),
                "DPYD", ("rs72549303",), "6744GA>CA", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*10", 1), HaplotypeCall("*2B", 1), HaplotypeCall("*9", 1)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("T", "T"), "DPYD", ("rs3918290",),
                "9213C>T", FullCallFilter.PASS, "9213C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TG"), "DPYD", ("rs72549303",),
                "REF_CALL", FullCallFilter.PASS, "6744GA>CA", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("C", "C"), "DPYD", ("rs1801159",),
                "293T>C", FullCallFilter.PASS, "293T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ambiguous_heterozygous_haplotype_without_a_clear_winner(self) -> None:
        """
        Ambiguous set of heterozygous haplotypes, with no haplotype that is simplest.
        Could be {*2A_HET,*10_HET}, {*5_HET,*9_HET} or {*3_HET,*2B_HET}.
        Uses v37 input.
        """
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("C", "T"),
                "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "293T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "T"), "DPYD", ("rs3918290",),
                "9213C>T", FullCallFilter.PASS, "9213C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "6744GA>CA", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "C"), "DPYD", ("rs1801159",),
                "293T>C", FullCallFilter.PASS, "293T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_unknown_variants_v37(self) -> None:
        """Variants that are completely unknown, including with unknown rs id. Uses v37 input"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("16", 39593405), "A"), ("A", "G"),
                "FAKE2", ("rs1949223",), "384C>T", SimpleCallFilter.PASS),  # unknown
            SimpleCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 2488242), "AC"), ("AC", "AG"),
                "DPYD", (".",), "9213CT>GT", SimpleCallFilter.PASS),  # unknown
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {
            "DPYD": set(), "FAKE": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)}, "FAKE2": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 2488242), "AC"), None,
                ("AC", "AG"), "DPYD", (".",),
                "9213CT>GT", FullCallFilter.PASS, "9213CT>GT?", FullCallFilter.UNKNOWN,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TG"), "DPYD", ("rs72549303",),
                "REF_CALL", FullCallFilter.NO_CALL, "6744GA>CA", FullCallFilter.INFERRED_PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 39593405), "A"), None,
                ("A", "G"), "FAKE2", ("rs1949223",),
                "384C>T", FullCallFilter.PASS, "384C>T?", FullCallFilter.UNKNOWN,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("C", "C"), "FAKE2", ("rs1212127",),
                "REF_CALL", FullCallFilter.NO_CALL, "1324T>C", FullCallFilter.INFERRED_PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "C"), "FAKE", ("rs1212125",),
                "1005T>C", FullCallFilter.PASS, "1005T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_unknown_variants_v38(self) -> None:
        """Variants that are completely unknown, including with unknown rs id. Uses v38 input"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr16", 39593405), "A"), ("A", "G"),
                "FAKE2", ("rs1949223",), "384C>T", SimpleCallFilter.PASS),  # unknown
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr5", 97450060), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 2488242), "AC"), ("AC", "AG"),
                "DPYD", (".",), "9213CT>GT", SimpleCallFilter.PASS),  # unknown
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {
            "DPYD": set(), "FAKE": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)}, "FAKE2": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                None, ReferenceSite(GeneCoordinate("chr1", 2488242), "AC"),
                ("AC", "AG"), "DPYD", (".",),
                "9213CT>GT?", FullCallFilter.UNKNOWN, "9213CT>GT", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.INFERRED_PASS, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                None, ReferenceSite(GeneCoordinate("chr16", 39593405), "A"),
                ("A", "G"), "FAKE2", ("rs1949223",),
                "384C>T?", FullCallFilter.UNKNOWN, "384C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("T", "T"), "FAKE2", ("rs1212127",),
                "1324C>T", FullCallFilter.INFERRED_PASS, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "C"), "FAKE", ("rs1212125",),
                "1005T>C", FullCallFilter.PASS, "1005T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_unknown_gene(self) -> None:
        """Variants that are of an unknown gene. Uses v37 input"""
        panel = self.__get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V37

        good_calls = {
            SimpleCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"), "FAKE",
                ("rs1212125",), "1005T>C", SimpleCallFilter.PASS),
        }
        bad_call = SimpleCall(
            ReferenceSite(GeneCoordinate("3", 18473423), "T"), ("T", "C"), ".",
            ("rs2492932",), "12T>C", SimpleCallFilter.PASS,
        )

        good_simple_call_data = SimpleCallData(frozenset(good_calls), reference_assembly)
        all_simple_call_data = SimpleCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser.create_pgx_analysis(good_simple_call_data, panel)  # no error
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(all_simple_call_data, panel)

    @unittest.skip("WIP")
    def test_incorrect_gene(self) -> None:
        """Variants that are assigned to an incorrect gene. Uses v37 input"""
        panel = self.__get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V37

        good_calls = {
            SimpleCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"), "FAKE",
                ("rs1212125",), "1005T>C", SimpleCallFilter.PASS),
        }
        bad_call = SimpleCall(
            ReferenceSite(GeneCoordinate("3", 18473423), "T"), ("T", "C"), "DPYD",
            ("rs2492932",), "12T>C", SimpleCallFilter.PASS,
        )

        good_simple_call_data = SimpleCallData(frozenset(good_calls), reference_assembly)
        all_simple_call_data = SimpleCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser.create_pgx_analysis(good_simple_call_data, panel)  # no error
        print(PgxAnalyser.create_pgx_analysis(all_simple_call_data, panel))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(all_simple_call_data, panel)

    def test_known_variants_with_incorrect_rs_id(self) -> None:
        """Known variant with one incorrect rs id. Uses v38 input"""
        panel = self.__get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V38

        good_calls = {
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "T"),
                "FAKE2", ("rs1212127",), "1324T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", SimpleCallFilter.PASS),
        }

        bad_call = SimpleCall(
            ReferenceSite(GeneCoordinate("chr5", 97450060), "T"), ("T", "C"),
            "FAKE", ("rs27384",), "1005T>C", SimpleCallFilter.PASS
        )

        good_call_data = SimpleCallData(frozenset(good_calls), reference_assembly)
        all_call_data = SimpleCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser.create_pgx_analysis(good_call_data, panel)
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(all_call_data, panel)

    def test_known_variant_with_incorrect_position(self) -> None:
        """Known variant with incorrect position. Uses v37 input."""
        panel = self.__get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V37

        good_calls = {
            SimpleCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("C", "T"),
                "FAKE2", ("rs1212127",), "1324C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", SimpleCallFilter.PASS),
        }
        bad_call = SimpleCall(
            ReferenceSite(GeneCoordinate("1", 6778543), "C"), ("C", "T"),
            "DPYD", ("rs3918290",), "35G>A", SimpleCallFilter.PASS
        )

        good_call_data = SimpleCallData(frozenset(good_calls), reference_assembly)
        all_call_data = SimpleCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser.create_pgx_analysis(good_call_data, panel)
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(all_call_data, panel)

    def test_known_variant_with_incorrect_chromosome(self) -> None:
        """Known variants with one incorrect chromosome. Uses v38 input"""
        panel = self.__get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V38

        good_calls = {
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "T"),
                "FAKE2", ("rs1212127",), "1324T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr5", 97450060), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744GA>CA", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", SimpleCallFilter.PASS),
        }
        bad_call = SimpleCall(
            ReferenceSite(GeneCoordinate("3", 97915614), "C"), ("C", "T"),
            "DPYD", ("rs3918290",), "35G>A", SimpleCallFilter.PASS
        )

        good_call_data = SimpleCallData(frozenset(good_calls), reference_assembly)
        all_call_data = SimpleCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser.create_pgx_analysis(good_call_data, panel)
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(all_call_data, panel)

    def test_known_variant_with_multiple_rs_ids_not_matching_panel(self) -> None:
        """Multiple rs ids when panel says there should be one. Uses v37 input."""
        panel = self.__get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V37
        good_calls = {
            SimpleCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", SimpleCallFilter.PASS),
        }
        bad_call = SimpleCall(
            ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("C", "T"),
            "FAKE2", ("rs1212127", "rs394832"), "1324C>T", SimpleCallFilter.PASS
        )

        good_call_data = SimpleCallData(frozenset(good_calls), reference_assembly)
        all_call_data = SimpleCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser.create_pgx_analysis(good_call_data, panel)
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(all_call_data, panel)

    def test_unresolved_haplotype_because_of_unexpected_base_v37(self) -> None:
        """No haplotype call because of unexpected base at known variant location, for v37 input."""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        reference_assembly = ReferenceAssembly.V37

        call_that_ref_seq_diff_is_ref_v38 = SimpleCall(
            ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TC", "TC"),
            "DPYD", ("rs72549303",), "6744CA>GA", SimpleCallFilter.PASS,
        )
        normal_call = SimpleCall(
            ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("T", "C"),
            "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS,
        )
        unexpected_base_call = SimpleCall(
            ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("T", "A"),
            "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS,
        )
        good_call_data = SimpleCallData(
            frozenset({normal_call, call_that_ref_seq_diff_is_ref_v38}),
            reference_assembly
        )
        unexpected_base_call_data = SimpleCallData(
            frozenset({unexpected_base_call, call_that_ref_seq_diff_is_ref_v38}),
            reference_assembly
        )
        good_pgx_analysis = PgxAnalyser.create_pgx_analysis(good_call_data, panel)
        good_gene_to_haplotype_calls_expected = {"DPYD": {HaplotypeCall("*2A", 1), HaplotypeCall("*1", 1)}}
        self.assertEqual(
            good_gene_to_haplotype_calls_expected, good_pgx_analysis.get_gene_to_haplotype_calls())

        unexpected_base_pgx_analysis = PgxAnalyser.create_pgx_analysis(unexpected_base_call_data, panel)

        unexpected_base_gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        unexpected_base_all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("T", "A"), "DPYD", ("rs3918290",),
                "9213C>T", FullCallFilter.PASS, "9213C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
        })
        unexpected_base_pgx_analysis_expected = PgxAnalysis(
            FullCallData(unexpected_base_all_full_calls_expected),
            unexpected_base_gene_to_haplotype_calls_expected)
        self.assertEqual(unexpected_base_pgx_analysis_expected, unexpected_base_pgx_analysis)

    def test_unresolved_haplotype_because_of_unexpected_base_v38(self) -> None:
        """No haplotype call because of unexpected base at known variant location, for v38 input."""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        reference_assembly = ReferenceAssembly.V38

        normal_call = SimpleCall(
            ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "C"),
            "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS,
        )
        unexpected_base_call = SimpleCall(
            ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "A"),
            "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS,
        )
        other_call = SimpleCall(
            ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("C", "C"),
            "DPYD", ("rs1801159",), "674A>G", SimpleCallFilter.PASS,
        )
        good_call_data = SimpleCallData(
            frozenset({normal_call, other_call}),
            reference_assembly
        )
        unexpected_base_call_data = SimpleCallData(
            frozenset({unexpected_base_call, other_call}),
            reference_assembly
        )
        good_pgx_analysis = PgxAnalyser.create_pgx_analysis(good_call_data, panel)
        good_gene_to_haplotype_calls_expected = {"DPYD": {HaplotypeCall("*2B", 1), HaplotypeCall("*5", 1)}}
        self.assertEqual(
            good_gene_to_haplotype_calls_expected, good_pgx_analysis.get_gene_to_haplotype_calls())

        unexpected_base_pgx_analysis = PgxAnalyser.create_pgx_analysis(unexpected_base_call_data, panel)

        unexpected_base_gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        unexpected_base_all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("T", "A"), "DPYD", ("rs3918290",),
                "9213C>T", FullCallFilter.PASS, "9213C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.INFERRED_PASS, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("C", "C"), "DPYD", ("rs1801159",),
                "674A>G", FullCallFilter.PASS, "674A>G", FullCallFilter.PASS,
            ),
        })
        unexpected_base_pgx_analysis_expected = PgxAnalysis(
            FullCallData(unexpected_base_all_full_calls_expected),
            unexpected_base_gene_to_haplotype_calls_expected)
        self.assertEqual(unexpected_base_pgx_analysis_expected, unexpected_base_pgx_analysis)

    def test_unresolved_haplotype_because_of_only_half_of_haplotype(self) -> None:
        """
        No haplotype call because of missing half of haplotype.
        One half is present twice, other once. Uses v37 data.
        """
        panel = self.__get_narrow_example_panel({"*2B"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9213C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "293T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("T", "T"), "DPYD", ("rs3918290",),
                "9213C>T", FullCallFilter.PASS, "9213C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "C"), "DPYD", ("rs1801159",),
                "293T>C", FullCallFilter.PASS, "293T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_unresolved_haplotype_because_mnv_covers_snv_starting_early(self) -> None:
        """
        Unresolved haplotype because MNV covers where SNV was expected, where MNV starts before this location.
        Uses v38 input.
        """
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450057), "GC"), ("CT", "CT"),
                "DPYD", (".",), "9212GC>CT", SimpleCallFilter.PASS),  # unknown
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "REF_CALL", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                None, ReferenceSite(GeneCoordinate("chr1", 97450057), "GC"),
                ("CT", "CT"), "DPYD", (".",),
                "9212GC>CT?", FullCallFilter.UNKNOWN, "9212GC>CT", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("C", "C"), "DPYD", ("rs1801159",),
                "293T>C", FullCallFilter.PASS, "293T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_unresolved_haplotype_because_mnv_covers_snv_starting_there(self) -> None:
        """
        Unresolved haplotype because MNV covers where SNV was expected, where MNV starts at this location.
        Uses v37 input.
        """
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "CG"), ("TC", "TC"),
                "DPYD", (".",), "9212CG>TC", SimpleCallFilter.PASS),  # unknown
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "CG"), None,
                ("TC", "TC"), "DPYD", (".",),
                "9212CG>TC", FullCallFilter.PASS, "9212CG>TC?", FullCallFilter.UNKNOWN,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("C", "C"), "DPYD", ("rs1801159",),
                "293T>C", FullCallFilter.PASS, "293T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ref_call_on_ref_seq_differences_v37(self) -> None:
        """Explicit ref calls wrt v37 at differences between v37 and v38"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("C", "C"),
                "FAKE2", ("rs1212127",), "REF_CALL", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TG"),
                "DPYD", ("rs72549303",), "REF_CALL", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*3", 2)}, "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": {HaplotypeCall("*4A", 2)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TG"), "DPYD", ("rs72549303",),
                "REF_CALL", FullCallFilter.PASS, "6744GA>CA", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("C", "C"), "FAKE2", ("rs1212127",),
                "REF_CALL", FullCallFilter.PASS, "1324T>C", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "T"), "FAKE", ("rs1212125",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ref_call_on_ref_seq_differences_v38(self) -> None:
        """Explicit ref calls wrt v38 at differences between v37 and v38"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("T", "T"),
                "FAKE2", ("rs1212127",), "REF_CALL", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "REF_CALL", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*1", 2)}, "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": {HaplotypeCall("*1", 2)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("T", "T"), "FAKE2", ("rs1212127",),
                "1324C>T", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "T"), "FAKE", ("rs1212125",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_homozygous_call_on_ref_seq_differences_v38(self) -> None:
        """Explicit ref calls wrt v38 at differences between v37 and v38"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "C"),
                "FAKE2", ("rs1212127",), "1324T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TG"),
                "DPYD", ("rs72549303",), "6744GA>CA", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*3", 2)}, "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": {HaplotypeCall("*4A", 2)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TG"), "DPYD", ("rs72549303",),
                "REF_CALL", FullCallFilter.PASS, "6744GA>CA", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("C", "C"), "FAKE2", ("rs1212127",),
                "REF_CALL", FullCallFilter.PASS, "1324T>C", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "T"), "FAKE", ("rs1212125",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_only_position_match_on_ref_seq_differences_v37(self) -> None:
        """At reference sequence differences: heterozygous between ref v37 and v38, and no rs_id provided"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("C", "T"),
                       "FAKE2", (".",), "1324C>T", SimpleCallFilter.PASS),
            SimpleCall(ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                       "DPYD", (".",), "6744CA>GA", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*3", 1), HaplotypeCall("*1", 1)},
            "FAKE": {HaplotypeCall("*1", 2)},
            "FAKE2": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)},
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "6744GA>CA", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("C", "T"), "FAKE2", ("rs1212127",),
                "1324C>T", FullCallFilter.PASS, "1324T>C", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "T"), "FAKE", ("rs1212125",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_only_position_match_on_ref_seq_differences_v38(self) -> None:
        """At reference sequence differences: heterozygous between ref v37 and v38, and no rs_id provided"""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "T"),
                       "FAKE2", (".",), "1324T>C", SimpleCallFilter.PASS),
            SimpleCall(ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TC"),
                       "DPYD", (".",), "6744GA>CA", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*3", 1), HaplotypeCall("*1", 1)},
            "FAKE": {HaplotypeCall("*1", 2)},
            "FAKE2": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)},
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "6744GA>CA", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("C", "T"), "FAKE2", ("rs1212127",),
                "1324C>T", FullCallFilter.PASS, "1324T>C", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "T"), "FAKE", ("rs1212125",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_single_different_allele_on_ref_seq_differences(self) -> None:
        """
        At reference sequence differences: single allele that is ref v37 or v38, other allele is neither.
        Uses v38 input.
        """
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "A"),
                "FAKE2", ("rs1212127",), "1324T>A", SimpleCallFilter.PASS),  # with ref v37
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("AC", "TC"),
                "DPYD", ("rs72549303",), "6744CT>GT;6744CT>GC", SimpleCallFilter.PASS),  # with ref v38
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {
            "DPYD": set(), "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": set()
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("AC", "TC"), "DPYD", ("rs72549303",),
                "6744CT>GT;6744CT>GC?", FullCallFilter.UNKNOWN, "6744CT>GT;6744CT>GC", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("C", "A"), "FAKE2", ("rs1212127",),
                "1324T>A?", FullCallFilter.UNKNOWN, "1324T>A", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "T"), "FAKE", ("rs1212125",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_double_different_allele_on_ref_seq_differences(self) -> None:
        """At reference sequence differences: both alleles not ref v37 or v38. Uses v37 input."""
        panel = self.__get_wide_example_panel()
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("A", "G"),
                "FAKE2", ("rs1212127",), "1324C>A;1324C>G", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("AC", "AG"),
                "DPYD", ("rs72549303",), "6744CT>GT;6744CT>GC", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {
            "DPYD": set(), "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": set()
        }
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("C", "C"), "DPYD", ("rs3918290",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("AC", "AG"), "DPYD", ("rs72549303",),
                "6744CT>GT;6744CT>GC", FullCallFilter.PASS, "6744CT>GT;6744CT>GC?", FullCallFilter.UNKNOWN,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("T", "T"), "DPYD", ("rs1801159",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"), ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"),
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ReferenceSite(GeneCoordinate("chr16", 97450060), "T"),
                ("A", "G"), "FAKE2", ("rs1212127",),
                "1324C>A;1324C>G", FullCallFilter.PASS, "1324C>A;1324C>G?", FullCallFilter.UNKNOWN,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ReferenceSite(GeneCoordinate("chr5", 97450060), "T"),
                ("T", "T"), "FAKE", ("rs1212125",),
                "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_overlapping_mnv_and_snv(self) -> None:
        """Overlapping SNV and unrecognized MNV do not cause errors. Uses v38 input."""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9212C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "CG"), ("CG", "TC"),
                "DPYD", (".",), "9212CG>TC", SimpleCallFilter.PASS),  # unknown
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "REF_CALL", SimpleCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("T", "T"), "DPYD", ("rs3918290",),
                "9212C>T", FullCallFilter.PASS, "9212C>T", FullCallFilter.PASS,
            ),
            FullCall(
                None, ReferenceSite(GeneCoordinate("chr1", 97450058), "CG"),
                ("CG", "TC"), "DPYD", (".",),
                "9212CG>TC?", FullCallFilter.UNKNOWN, "9212CG>TC", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TC", "TC"), "DPYD", ("rs72549303",),
                "6744CA>GA", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("C", "C"), "DPYD", ("rs1801159",),
                "293T>C", FullCallFilter.PASS, "293T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_overlapping_mnv_and_snv_at_ref_seq_difference(self) -> None:
        """Overlapping SNV and unrecognized MNV at ref seq difference do not cause error. Uses v37 input."""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        simple_call_data = SimpleCallData(frozenset({
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9212C>T", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "AC"),
                "DPYD", ("rs72549303",), "6744CA>GT", SimpleCallFilter.PASS),
            SimpleCall(
                ReferenceSite(GeneCoordinate("1", 97915622), "G"), ("C", "C"),
                "DPYD", (".",), "6744C>G", SimpleCallFilter.PASS),  # unknown
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser.create_pgx_analysis(simple_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("T", "T"), "DPYD", ("rs3918290",),
                "9212C>T", FullCallFilter.PASS, "9212C>T", FullCallFilter.PASS,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"),
                ("TG", "AC"), "DPYD", ("rs72549303",),
                "6744CA>GT", FullCallFilter.PASS, '6744CA>GT?', FullCallFilter.UNKNOWN,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915622), "G"), None,
                ("C", "C"), "DPYD", (".",),
                "6744C>G", FullCallFilter.PASS, "6744C>G?", FullCallFilter.UNKNOWN,
            ),
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ReferenceSite(GeneCoordinate("chr1", 97515839), "T"),
                ("C", "C"), "DPYD", ("rs1801159",),
                "293T>C", FullCallFilter.PASS, "293T>C", FullCallFilter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(FullCallData(all_full_calls_expected), gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
