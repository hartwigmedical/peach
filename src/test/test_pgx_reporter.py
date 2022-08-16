import unittest
from typing import Dict, Set

from util.filter import DualCallFilter
from util.gene_coordinate import GeneCoordinate
from util.reference_assembly import ReferenceAssembly
from util.reference_site import ReferenceSite
from calls.haplotype_call import HaplotypeCall
from calls.dual_call import DualCall, DualCallData
from analysis.pgx_analysis import PgxAnalysis
from pgx_reporter import GenotypeReporter, HaplotypeReporter
from test.util_for_test import get_wide_example_panel, get_empty_panel


class TestPgxReporter(unittest.TestCase):
    def test_genotype_reporter_empty(self) -> None:
        pgx_analysis = PgxAnalysis(DualCallData(frozenset()), {})
        panel_id = "Panel_v0.2"
        version = "V1"
        result = GenotypeReporter().get_calls_tsv_text(pgx_analysis, panel_id, version, ReferenceAssembly.V37)

        result_expected = (
            "gene\tchromosome\tposition_v37\tposition_v38\tref_v37\tref_v38\t"
            "allele1\tallele2\trsid\tvariant_annotation_v37\tfilter_v37\tvariant_annotation_v38\tfilter_v38\tpanel_version\trepo_version\tchromosome_v38\n"
        )
        self.assertEqual(result_expected, result)

    def test_genotype_reporter_non_empty_input_v37(self) -> None:
        all_dual_calls = frozenset({
            DualCall(
                ReferenceSite(GeneCoordinate("1", 1), "C"), None, ("C", "CAG"),
                "DPYD", ("rs664",), "1A>C;1A>G", DualCallFilter.PASS, "1A>C;1A>G?", DualCallFilter.UNKNOWN,
            ),
            DualCall(
                ReferenceSite(GeneCoordinate("1", 5), "A"), ReferenceSite(GeneCoordinate("chr1", 25), "A"), ("G", "C"),
                "DPYD", tuple(), None, DualCallFilter.PASS, None, DualCallFilter.PASS,
            ),
            DualCall(
                ReferenceSite(GeneCoordinate("1", 15), "C"), None, ("C", "CAG"),
                "DPYD", ("rs536",), "35A>C;35A>G", DualCallFilter.PASS, "35A>C;35A>G?", DualCallFilter.UNKNOWN,
            ),
            DualCall(
                ReferenceSite(GeneCoordinate("X", 15), "TT"), ReferenceSite(GeneCoordinate("chrX", 40), "AA"), ("TT", "TT"),
                "GENE", ("rs23",), "REF_CALL", DualCallFilter.NO_CALL, "627AA>TT", DualCallFilter.INFERRED_PASS,
            ),
            DualCall(
                ReferenceSite(GeneCoordinate("2", 154663), "T"), ReferenceSite(GeneCoordinate("chr2", 40565464), "T"), ("T", "T"),
                "BRAF", ("rs154", "rs8839"), "REF_CALL", DualCallFilter.NO_CALL, "REF_CALL", DualCallFilter.NO_CALL,
            ),
            DualCall(
                ReferenceSite(GeneCoordinate("15", 24113), "A"), ReferenceSite(GeneCoordinate("chr15", 684633), "T"), ("T", "T"),
                ".", ("rs462", "rs9820", "rs536"), "29482A>T", DualCallFilter.PASS, "REF_CALL", DualCallFilter.PASS,
            ),
        })
        pgx_analysis = PgxAnalysis(DualCallData(all_dual_calls), {})
        panel_id = "Panel_v0.2"
        version = "V1"
        result = GenotypeReporter().get_calls_tsv_text(pgx_analysis, panel_id, version, ReferenceAssembly.V37)

        result_expected = (
            "gene\tchromosome\tposition_v37\tposition_v38\tref_v37\tref_v38\tallele1\tallele2\t"
                "rsid\tvariant_annotation_v37\tfilter_v37\tvariant_annotation_v38\tfilter_v38\tpanel_version\trepo_version\tchromosome_v38\n"
            "DPYD\t1\t1\tUNKNOWN\tC\tUNKNOWN\tC\tCAG\trs664\t1A>C;1A>G\tPASS\t1A>C;1A>G?\tUNKNOWN\tPanel_v0.2\tV1\tUNKNOWN\n"
            "DPYD\t1\t5\t25\tA\tA\tC\tG\t.\tUNKNOWN\tPASS\tUNKNOWN\tPASS\tPanel_v0.2\tV1\tchr1\n"
            "DPYD\t1\t15\tUNKNOWN\tC\tUNKNOWN\tC\tCAG\trs536\t35A>C;35A>G\tPASS\t35A>C;35A>G?\tUNKNOWN\tPanel_v0.2\tV1\tUNKNOWN\n"
            "BRAF\t2\t154663\t40565464\tT\tT\tT\tT\trs154;rs8839\tREF_CALL\tNO_CALL\tREF_CALL\tNO_CALL\tPanel_v0.2\tV1\tchr2\n"
            ".\t15\t24113\t684633\tA\tT\tT\tT\trs462;rs9820;rs536\t29482A>T\tPASS\tREF_CALL\tPASS\tPanel_v0.2\tV1\tchr15\n"
            "GENE\tX\t15\t40\tTT\tAA\tTT\tTT\trs23\tREF_CALL\tNO_CALL\t627AA>TT\tINFERRED_PASS\tPanel_v0.2\tV1\tchrX\n"
        )
        self.assertEqual(result_expected, result)

    def test_genotype_reporter_non_empty_input_v38(self) -> None:
        all_dual_calls = frozenset({
            DualCall(
                None, ReferenceSite(GeneCoordinate("chr1", 15), "C"), ("C", "CAG"),
                "DPYD", ("rs353",), "3A>C;3A>G?", DualCallFilter.UNKNOWN, "3A>C;3A>G", DualCallFilter.PASS,
            ),
            DualCall(
                ReferenceSite(GeneCoordinate("1", 5), "A"), ReferenceSite(GeneCoordinate("chr1", 25), "A"), ("G", "C"),
                "DPYD", tuple(), None, DualCallFilter.PASS, None, DualCallFilter.PASS,
            ),
            DualCall(
                None, ReferenceSite(GeneCoordinate("chr1", 35), "C"), ("C", "CAG"),
                "DPYD", ("rs536",), "35A>C;35A>G?", DualCallFilter.UNKNOWN, "35A>C;35A>G", DualCallFilter.PASS,
            ),
            DualCall(
                ReferenceSite(GeneCoordinate("X", 15), "AA"), ReferenceSite(GeneCoordinate("chrX", 40), "TT"), ("TT", "TT"),
                "GENE", ("rs23",), "627AA>TT", DualCallFilter.INFERRED_PASS, "REF_CALL", DualCallFilter.NO_CALL,
            ),
            DualCall(
                ReferenceSite(GeneCoordinate("2", 154663), "T"), ReferenceSite(GeneCoordinate("chr2", 40565464), "T"), ("T", "T"),
                "BRAF", ("rs154", "rs8839"), "REF_CALL", DualCallFilter.NO_CALL, "REF_CALL", DualCallFilter.NO_CALL,
            ),
            DualCall(
                ReferenceSite(GeneCoordinate("15", 24113), "T"), ReferenceSite(GeneCoordinate("chr15", 684633), "A"), ("T", "T"),
                ".", ("rs462", "rs9820", "rs536"), "REF_CALL", DualCallFilter.PASS, "29482A>T", DualCallFilter.PASS,
            ),
        })
        pgx_analysis = PgxAnalysis(DualCallData(all_dual_calls), {})
        panel_id = "Panel_v0.2"
        version = "V1"
        result = GenotypeReporter().get_calls_tsv_text(pgx_analysis, panel_id, version, ReferenceAssembly.V38)

        result_expected = (
            "gene\tchromosome\tposition_v37\tposition_v38\tref_v37\tref_v38\tallele1\tallele2\t"
                "rsid\tvariant_annotation_v37\tfilter_v37\tvariant_annotation_v38\tfilter_v38\tpanel_version\trepo_version\tchromosome_v38\n"
            "DPYD\tUNKNOWN\tUNKNOWN\t15\tUNKNOWN\tC\tC\tCAG\trs353\t3A>C;3A>G?\tUNKNOWN\t3A>C;3A>G\tPASS\tPanel_v0.2\tV1\tchr1\n"
            "DPYD\t1\t5\t25\tA\tA\tC\tG\t.\tUNKNOWN\tPASS\tUNKNOWN\tPASS\tPanel_v0.2\tV1\tchr1\n"
            "DPYD\tUNKNOWN\tUNKNOWN\t35\tUNKNOWN\tC\tC\tCAG\trs536\t35A>C;35A>G?\tUNKNOWN\t35A>C;35A>G\tPASS\tPanel_v0.2\tV1\tchr1\n"
            "BRAF\t2\t154663\t40565464\tT\tT\tT\tT\trs154;rs8839\tREF_CALL\tNO_CALL\tREF_CALL\tNO_CALL\tPanel_v0.2\tV1\tchr2\n"
            ".\t15\t24113\t684633\tT\tA\tT\tT\trs462;rs9820;rs536\tREF_CALL\tPASS\t29482A>T\tPASS\tPanel_v0.2\tV1\tchr15\n"
            "GENE\tX\t15\t40\tAA\tTT\tTT\tTT\trs23\t627AA>TT\tINFERRED_PASS\tREF_CALL\tNO_CALL\tPanel_v0.2\tV1\tchrX\n"
        )
        self.assertEqual(result_expected, result)

    def test_haplotype_reporter_empty(self) -> None:
        pgx_analysis = PgxAnalysis(DualCallData(frozenset()), {})
        panel = get_empty_panel()
        version = "V1"
        result = HaplotypeReporter().get_genotype_tsv_text(pgx_analysis, panel, version)

        result_expected = (
            "gene\thaplotype\tfunction\tlinked_drugs\turl_prescription_info\tpanel_version\trepo_version\thaplotype_only\tzygosity_only\n"
        )
        self.assertEqual(result_expected, result)

    def test_haplotype_reporter_non_empty(self) -> None:
        gene_to_haplotype_calls: Dict[str, Set[HaplotypeCall]] = {
            "DPYD": {HaplotypeCall("*2B", 2), HaplotypeCall("*3", 1), HaplotypeCall("*2A", 1)},
            "FAKE": set(),
            "FAKE2": {HaplotypeCall("*1", 1), HaplotypeCall("*4A", 1)},
        }
        pgx_analysis = PgxAnalysis(DualCallData(frozenset()), gene_to_haplotype_calls)
        panel = get_wide_example_panel(include_transcript_ids=True)
        version = "V1"
        result = HaplotypeReporter().get_genotype_tsv_text(pgx_analysis, panel, version)

        result_expected = (
            "gene\thaplotype\tfunction\tlinked_drugs\turl_prescription_info\tpanel_version\trepo_version\thaplotype_only\tzygosity_only\n"
            "DPYD\t*2A_HET\tNo Function\t5-Fluorouracil;Capecitabine\thttps://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963\tWideTestPanel_v1.0\tV1\t*2A\tHET\n"
            "DPYD\t*2B_HOM\tNo Function\t5-Fluorouracil;Capecitabine\thttps://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963\tWideTestPanel_v1.0\tV1\t*2B\tHOM\n"
            "DPYD\t*3_HET\tNormal Function\t5-Fluorouracil;Capecitabine\thttps://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963\tWideTestPanel_v1.0\tV1\t*3\tHET\n"
            "FAKE\tUnresolved Haplotype\tUnknown Function\tAspirin\thttps://www.pharmgkb.org/some_other_url\tWideTestPanel_v1.0\tV1\tUnresolved Haplotype\tN/A\n"
            "FAKE2\t*1_HET\tNormal Function\tAspirin\thttps://www.pharmgkb.org/some_other_url\tWideTestPanel_v1.0\tV1\t*1\tHET\n"
            "FAKE2\t*4A_HET\tReduced Function\tAspirin\thttps://www.pharmgkb.org/some_other_url\tWideTestPanel_v1.0\tV1\t*4A\tHET\n"
        )
        self.assertEqual(result_expected, result)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
