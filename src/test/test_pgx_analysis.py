import unittest
from typing import Dict, Set

from base.filter import FullCallFilter, VcfCallFilter
from base.gene_coordinate import GeneCoordinate
from base.reference_assembly import ReferenceAssembly
from base.reference_site import ReferenceSite
from call_data import VcfCallData, VcfCall, HaplotypeCall, FullCall, FullCallData
from analysis.pgx_analysis import PgxAnalyser, PgxAnalysis
from test.util_for_test import get_wide_example_panel, get_narrow_example_panel


class TestPgxAnalysis(unittest.TestCase):
    def test_empty_v37(self) -> None:
        """No variants wrt v37"""
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset(), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset(), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("T", "T"),
                "FAKE2", ("rs1212127",), "1324C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TC", "TC"),
                "DPYD", tuple(), "6744CA>GA", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("C", "C"),
                "DPYD", ("rs3918290",), "REF_CALL", VcfCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("C", "T"),
                "FAKE2", ("rs1212127",), "1324C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("C", "T"),
                "DPYD", ("rs3918290",), "35G>A", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", VcfCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "T"),
                "FAKE2", ("rs1212127",), "1324T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr5", 97450060), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744GA>CA", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("C", "T"),
                "DPYD", ("rs3918290",), "35G>A", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", VcfCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_narrow_example_panel({"*2A", "*3", "*7"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("C", "T"),
                "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97912838), "A"), ("AGT", "AGT"),
                "DPYD", ("rs2938101",), "293A>AGT", VcfCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_narrow_example_panel({"*2A", "*5", "*2B"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "REF_CALL", VcfCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_narrow_example_panel({"*2A", "*5", "*2B"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("C", "T"),
                "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "293T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", VcfCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_narrow_example_panel({"*2A", "*5", "*2B"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "293T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "REF_CALL", VcfCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10", "*7"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "293T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97912838), "A"), ("AGT", "AGT"),
                "DPYD", ("rs2938101",), "301A>AGT", VcfCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TG"),
                "DPYD", ("rs72549303",), "6744GA>CA", VcfCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("C", "T"),
                "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "293T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", VcfCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        """
        Variants that are completely unknown, including with unknown rs id. Uses v37 input

        This should never happen because these calls should not be extracted from the VCF in the first place.
        """
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("16", 39593405), "A"), ("A", "G"),
                "FAKE2", ("rs1949223",), "384C>T", VcfCallFilter.PASS),  # unknown
            VcfCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 2488242), "AC"), ("AC", "AG"),
                "DPYD", tuple(), "9213CT>GT", VcfCallFilter.PASS),  # unknown
        }), ReferenceAssembly.V37)

        with self.assertRaises(ValueError):
            PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

    def test_unknown_variants_v38(self) -> None:
        """
        Variants that are completely unknown, including with unknown rs id. Uses v38 input.

        This should never happen because these calls should not be extracted from the VCF in the first place.
        """
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("chr16", 39593405), "A"), ("A", "G"),
                "FAKE2", ("rs1949223",), "384C>T", VcfCallFilter.PASS),  # unknown
            VcfCall(
                ReferenceSite(GeneCoordinate("chr5", 97450060), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 2488242), "AC"), ("AC", "AG"),
                "DPYD", tuple(), "9213CT>GT", VcfCallFilter.PASS),  # unknown
        }), ReferenceAssembly.V38)
        with self.assertRaises(ValueError):
            PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

    def test_unknown_gene(self) -> None:
        """Variants that are of an unknown gene. Uses v37 input"""
        panel = get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V37

        good_calls = {
            VcfCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"), "FAKE",
                ("rs1212125",), "1005T>C", VcfCallFilter.PASS),
        }
        bad_call = VcfCall(
            ReferenceSite(GeneCoordinate("3", 18473423), "T"), ("T", "C"), ".",
            ("rs2492932",), "12T>C", VcfCallFilter.PASS,
        )

        good_vcf_call_data = VcfCallData(frozenset(good_calls), reference_assembly)
        all_vcf_call_data = VcfCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser().create_pgx_analysis(good_vcf_call_data, panel)  # no error
        with self.assertRaises(ValueError):
            PgxAnalyser().create_pgx_analysis(all_vcf_call_data, panel)

    def test_incorrect_gene(self) -> None:
        """Variants that are assigned to an incorrect gene. Uses v37 input"""
        panel = get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V37

        good_calls = {
            VcfCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"), "FAKE",
                ("rs1212125",), "1005T>C", VcfCallFilter.PASS),
        }
        bad_call = VcfCall(
            ReferenceSite(GeneCoordinate("3", 18473423), "T"), ("T", "C"), "DPYD",
            ("rs2492932",), "12T>C", VcfCallFilter.PASS,
        )

        good_vcf_call_data = VcfCallData(frozenset(good_calls), reference_assembly)
        all_vcf_call_data = VcfCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser().create_pgx_analysis(good_vcf_call_data, panel)  # no error
        with self.assertRaises(ValueError):
            PgxAnalyser().create_pgx_analysis(all_vcf_call_data, panel)

    def test_known_variants_with_incorrect_rs_id(self) -> None:
        """Known variant with one incorrect rs id. Uses v38 input"""
        panel = get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V38

        good_calls = {
            VcfCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "T"),
                "FAKE2", ("rs1212127",), "1324T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", VcfCallFilter.PASS),
        }

        bad_call = VcfCall(
            ReferenceSite(GeneCoordinate("chr5", 97450060), "T"), ("T", "C"),
            "FAKE", ("rs27384",), "1005T>C", VcfCallFilter.PASS
        )

        good_call_data = VcfCallData(frozenset(good_calls), reference_assembly)
        all_call_data = VcfCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser().create_pgx_analysis(good_call_data, panel)
        with self.assertRaises(ValueError):
            PgxAnalyser().create_pgx_analysis(all_call_data, panel)

    def test_known_variant_with_incorrect_position(self) -> None:
        """Known variant with incorrect position. Uses v37 input."""
        panel = get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V37

        good_calls = {
            VcfCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("C", "T"),
                "FAKE2", ("rs1212127",), "1324C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", VcfCallFilter.PASS),
        }
        bad_call = VcfCall(
            ReferenceSite(GeneCoordinate("1", 6778543), "C"), ("C", "T"),
            "DPYD", ("rs3918290",), "35G>A", VcfCallFilter.PASS
        )

        good_call_data = VcfCallData(frozenset(good_calls), reference_assembly)
        all_call_data = VcfCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser().create_pgx_analysis(good_call_data, panel)
        with self.assertRaises(ValueError):
            PgxAnalyser().create_pgx_analysis(all_call_data, panel)

    def test_known_variant_with_incorrect_chromosome(self) -> None:
        """Known variants with one incorrect chromosome. Uses v38 input"""
        panel = get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V38

        good_calls = {
            VcfCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "T"),
                "FAKE2", ("rs1212127",), "1324T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr5", 97450060), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744GA>CA", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", VcfCallFilter.PASS),
        }
        bad_call = VcfCall(
            ReferenceSite(GeneCoordinate("3", 97915614), "C"), ("C", "T"),
            "DPYD", ("rs3918290",), "35G>A", VcfCallFilter.PASS
        )

        good_call_data = VcfCallData(frozenset(good_calls), reference_assembly)
        all_call_data = VcfCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser().create_pgx_analysis(good_call_data, panel)
        with self.assertRaises(ValueError):
            PgxAnalyser().create_pgx_analysis(all_call_data, panel)

    def test_known_variant_with_multiple_rs_ids_not_matching_panel(self) -> None:
        """Multiple rs ids when panel says there should be one. Uses v37 input."""
        panel = get_wide_example_panel()
        reference_assembly = ReferenceAssembly.V37
        good_calls = {
            VcfCall(
                ReferenceSite(GeneCoordinate("5", 97915617), "T"), ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", VcfCallFilter.PASS),
        }
        bad_call = VcfCall(
            ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("C", "T"),
            "FAKE2", ("rs1212127", "rs394832"), "1324C>T", VcfCallFilter.PASS
        )

        good_call_data = VcfCallData(frozenset(good_calls), reference_assembly)
        all_call_data = VcfCallData(frozenset(good_calls.union({bad_call})), reference_assembly)

        PgxAnalyser().create_pgx_analysis(good_call_data, panel)
        with self.assertRaises(ValueError):
            PgxAnalyser().create_pgx_analysis(all_call_data, panel)

    def test_unresolved_haplotype_because_of_unexpected_base_v37(self) -> None:
        """No haplotype call because of unexpected base at known variant location, for v37 input."""
        panel = get_narrow_example_panel({"*2A", "*5", "*2B"})
        reference_assembly = ReferenceAssembly.V37

        call_that_ref_seq_diff_is_ref_v38 = VcfCall(
            ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TC", "TC"),
            "DPYD", ("rs72549303",), "6744CA>GA", VcfCallFilter.PASS,
        )
        normal_call = VcfCall(
            ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("T", "C"),
            "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS,
        )
        unexpected_base_call = VcfCall(
            ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("T", "A"),
            "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS,
        )
        good_call_data = VcfCallData(
            frozenset({normal_call, call_that_ref_seq_diff_is_ref_v38}),
            reference_assembly
        )
        unexpected_base_call_data = VcfCallData(
            frozenset({unexpected_base_call, call_that_ref_seq_diff_is_ref_v38}),
            reference_assembly
        )
        good_pgx_analysis = PgxAnalyser().create_pgx_analysis(good_call_data, panel)
        good_gene_to_haplotype_calls_expected = {"DPYD": {HaplotypeCall("*2A", 1), HaplotypeCall("*1", 1)}}
        self.assertEqual(
            good_gene_to_haplotype_calls_expected, good_pgx_analysis.get_gene_to_haplotype_calls())

        unexpected_base_pgx_analysis = PgxAnalyser().create_pgx_analysis(unexpected_base_call_data, panel)

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
        panel = get_narrow_example_panel({"*2A", "*5", "*2B"})
        reference_assembly = ReferenceAssembly.V38

        normal_call = VcfCall(
            ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "C"),
            "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS,
        )
        unexpected_base_call = VcfCall(
            ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "A"),
            "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS,
        )
        other_call = VcfCall(
            ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("C", "C"),
            "DPYD", ("rs1801159",), "674A>G", VcfCallFilter.PASS,
        )
        good_call_data = VcfCallData(
            frozenset({normal_call, other_call}),
            reference_assembly
        )
        unexpected_base_call_data = VcfCallData(
            frozenset({unexpected_base_call, other_call}),
            reference_assembly
        )
        good_pgx_analysis = PgxAnalyser().create_pgx_analysis(good_call_data, panel)
        good_gene_to_haplotype_calls_expected = {"DPYD": {HaplotypeCall("*2B", 1), HaplotypeCall("*5", 1)}}
        self.assertEqual(
            good_gene_to_haplotype_calls_expected, good_pgx_analysis.get_gene_to_haplotype_calls())

        unexpected_base_pgx_analysis = PgxAnalyser().create_pgx_analysis(unexpected_base_call_data, panel)

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
        panel = get_narrow_example_panel({"*2B"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9213C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("T", "C"),
                "DPYD", ("rs1801159",), "293T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", VcfCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_narrow_example_panel({"*2A", "*5", "*2B"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450057), "GC"), ("CT", "CT"),
                "DPYD", tuple(), "9212GC>CT", VcfCallFilter.PASS),  # unknown
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "REF_CALL", VcfCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                None, ReferenceSite(GeneCoordinate("chr1", 97450057), "GC"),
                ("CT", "CT"), "DPYD", tuple(),
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
        panel = get_narrow_example_panel({"*2A", "*5", "*2B"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "CG"), ("TC", "TC"),
                "DPYD", tuple(), "9212CG>TC", VcfCallFilter.PASS),  # unknown
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", VcfCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "CG"), None,
                ("TC", "TC"), "DPYD", tuple(),
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
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("C", "C"),
                "FAKE2", ("rs1212127",), "REF_CALL", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TG"),
                "DPYD", ("rs72549303",), "REF_CALL", VcfCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("T", "T"),
                "FAKE2", ("rs1212127",), "REF_CALL", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "REF_CALL", VcfCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "C"),
                "FAKE2", ("rs1212127",), "1324T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TG"),
                "DPYD", ("rs72549303",), "6744GA>CA", VcfCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("C", "T"),
                       "FAKE2", tuple(), "1324C>T", VcfCallFilter.PASS),
            VcfCall(ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "TC"),
                       "DPYD", tuple(), "6744CA>GA", VcfCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "T"),
                       "FAKE2", tuple(), "1324T>C", VcfCallFilter.PASS),
            VcfCall(ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TG", "TC"),
                       "DPYD", tuple(), "6744GA>CA", VcfCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), ("C", "A"),
                "FAKE2", ("rs1212127",), "1324T>A", VcfCallFilter.PASS),  # with ref v37
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("AC", "TC"),
                "DPYD", ("rs72549303",), "6744CT>GT;6744CT>GC", VcfCallFilter.PASS),  # with ref v38
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_wide_example_panel()
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("16", 97915617), "C"), ("A", "G"),
                "FAKE2", ("rs1212127",), "1324C>A;1324C>G", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("AC", "AG"),
                "DPYD", ("rs72549303",), "6744CT>GT;6744CT>GC", VcfCallFilter.PASS),
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
        panel = get_narrow_example_panel({"*2A", "*5", "*2B"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9212C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450058), "CG"), ("CG", "TC"),
                "DPYD", tuple(), "9212CG>TC", VcfCallFilter.PASS),  # unknown
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), ("TC", "TC"),
                "DPYD", ("rs72549303",), "REF_CALL", VcfCallFilter.PASS),
        }), ReferenceAssembly.V38)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ReferenceSite(GeneCoordinate("chr1", 97450058), "C"),
                ("T", "T"), "DPYD", ("rs3918290",),
                "9212C>T", FullCallFilter.PASS, "9212C>T", FullCallFilter.PASS,
            ),
            FullCall(
                None, ReferenceSite(GeneCoordinate("chr1", 97450058), "CG"),
                ("CG", "TC"), "DPYD", tuple(),
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
        panel = get_narrow_example_panel({"*2A", "*5", "*2B"})
        vcf_call_data = VcfCallData(frozenset({
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915614), "C"), ("T", "T"),
                "DPYD", ("rs3918290",), "9212C>T", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97981395), "T"), ("C", "C"),
                "DPYD", ("rs1801159",), "293T>C", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915621), "TG"), ("TG", "AC"),
                "DPYD", ("rs72549303",), "6744CA>GT", VcfCallFilter.PASS),
            VcfCall(
                ReferenceSite(GeneCoordinate("1", 97915622), "G"), ("C", "C"),
                "DPYD", tuple(), "6744C>G", VcfCallFilter.PASS),  # unknown
        }), ReferenceAssembly.V37)
        pgx_analysis = PgxAnalyser().create_pgx_analysis(vcf_call_data, panel)

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
                ("C", "C"), "DPYD", tuple(),
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
