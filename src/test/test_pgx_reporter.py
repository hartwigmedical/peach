import unittest
from typing import Dict, Set

from base.filter import FullCallFilter
from base.gene_coordinate import GeneCoordinate
from base.reference_assembly import ReferenceAssembly
from call_data import FullCall, HaplotypeCall, FullCallData
from config.annotation import Annotation
from config.drug_info import DrugInfo
from config.gene_info import GeneInfo
from config.haplotype import Haplotype
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant
from pgx_analysis import PgxAnalysis
from pgx_reporter import GenotypeReporter, HaplotypeReporter


class TestPgxReporter(unittest.TestCase):
    @classmethod
    def __get_example_panel(cls) -> Panel:
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
            RsIdInfo("rs3918290", "C", "C", GeneCoordinate("1", 97915614), GeneCoordinate("chr1", 97450058)),
            RsIdInfo("rs72549309", "GATGA", "GATGA", GeneCoordinate("1", 98205966), GeneCoordinate("chr1", 97740410)),
            RsIdInfo("rs1801159", "T", "T", GeneCoordinate("1", 97981395), GeneCoordinate("chr1", 97515839)),
            RsIdInfo("rs72549303", "TG", "TC", GeneCoordinate("1", 97915621), GeneCoordinate("chr1", 97450065)),
        })
        dpyd_drugs = frozenset({
            DrugInfo("5-Fluorouracil", "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"),
            DrugInfo("Capecitabine", "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"),
        })
        dpyd_rs_id_to_difference_annotations = {
            "rs72549303": Annotation("6744CA>GA", "6744GA>CA"),
        }

        fake_haplotypes = frozenset({
            Haplotype("*4A", "Reduced Function", frozenset({fake_variant})),
        })
        fake_rs_id_infos = frozenset({
            RsIdInfo("rs1212125", "T", "T", GeneCoordinate("5", 97915617), GeneCoordinate("chr5", 97450060)),
        })
        fake_drugs = frozenset({
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
        })
        fake_rs_id_to_difference_annotations: Dict[str, Annotation] = {}

        fake2_haplotypes = frozenset({
            Haplotype("*4A", "Reduced Function", frozenset({fake2_variant})),
        })
        fake2_rs_id_infos = frozenset({
            RsIdInfo("rs1212127", "C", "T", GeneCoordinate("16", 97915617), GeneCoordinate("chr16", 97450060)),
        })
        fake2_drugs = frozenset({
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
        })
        fake2_rs_id_to_difference_annotations = {"rs1212127": Annotation("1324C>T", "1324T>C")}

        gene_infos = frozenset({
            GeneInfo("DPYD", "*1", dpyd_haplotypes, dpyd_rs_id_infos,
                     dpyd_drugs, dpyd_rs_id_to_difference_annotations),
            GeneInfo("FAKE", "*1", fake_haplotypes, fake_rs_id_infos,
                     fake_drugs, fake_rs_id_to_difference_annotations),
            GeneInfo("FAKE2", "*1", fake2_haplotypes, fake2_rs_id_infos,
                     fake2_drugs, fake2_rs_id_to_difference_annotations),
        })
        name = "Panel"
        version = "0.2"
        return Panel(name, version, gene_infos)

    def test_genotype_reporter_empty(self) -> None:
        pgx_analysis = PgxAnalysis(FullCallData(frozenset()), {})
        panel_id = "Panel_v0.2"
        version = "V1"
        result = GenotypeReporter.get_calls_tsv_text(pgx_analysis, panel_id, version, ReferenceAssembly.V37)

        result_expected = (
            "gene\tchromosome_v37\tposition_v37\tchromosome_v38\tposition_v38\tref_v37\tref_v38\t"
            "allele1\tallele2\trsid\tvariant_annotation_v37\tfilter_v37\tvariant_annotation_v38\tfilter_v38\tpanel_version\trepo_version\n"
        )
        self.assertEqual(result_expected, result)

    def test_genotype_reporter_non_empty_input_v37(self) -> None:
        all_full_calls = frozenset({
            FullCall(GeneCoordinate("1", 1), "C", None, None, ("C", "CAG"),
                     "DPYD", ("rs664",), "1A>C;1A>G", FullCallFilter.PASS, "1A>C;1A>G?", FullCallFilter.UNKNOWN),
            FullCall(GeneCoordinate("1", 5), "A", GeneCoordinate("chr1", 25), "A", ("G", "C"),
                     "DPYD", (".",), "25A>C;25A>G", FullCallFilter.PASS, "25A>C;25A>G", FullCallFilter.PASS),
            FullCall(GeneCoordinate("1", 15), "C", None, None, ("C", "CAG"),
                     "DPYD", ("rs536",), "35A>C;35A>G", FullCallFilter.PASS, "35A>C;35A>G?", FullCallFilter.UNKNOWN),
            FullCall(GeneCoordinate("X", 15), "TT", GeneCoordinate("chrX", 40), "AA", ("TT", "TT"),
                     "GENE", ("rs23",), "REF_CALL", FullCallFilter.NO_CALL, "627AA>TT", FullCallFilter.INFERRED_PASS),
            FullCall(GeneCoordinate("2", 154663), "T", GeneCoordinate("chr2", 40565464), "T", ("T", "T"),
                     "BRAF", ("rs154", "rs8839"), "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL),
            FullCall(GeneCoordinate("15", 24113), "A", GeneCoordinate("chr15", 684633), "T", ("T", "T"),
                     ".", ("rs462", "rs9820", "rs536"), "29482A>T", FullCallFilter.PASS, "REF_CALL", FullCallFilter.PASS),
        })
        pgx_analysis = PgxAnalysis(FullCallData(all_full_calls), {})
        panel_id = "Panel_v0.2"
        version = "V1"
        result = GenotypeReporter.get_calls_tsv_text(pgx_analysis, panel_id, version, ReferenceAssembly.V37)

        result_expected = (
            "gene\tchromosome_v37\tposition_v37\tchromosome_v38\tposition_v38\tref_v37\tref_v38\tallele1\tallele2\t"
                "rsid\tvariant_annotation_v37\tfilter_v37\tvariant_annotation_v38\tfilter_v38\tpanel_version\trepo_version\n"
            "DPYD\t1\t1\tUNKNOWN\tUNKNOWN\tC\tUNKNOWN\tC\tCAG\trs664\t1A>C;1A>G\tPASS\t1A>C;1A>G?\tUNKNOWN\tPanel_v0.2\tV1\n"
            "DPYD\t1\t5\tchr1\t25\tA\tA\tC\tG\t.\t25A>C;25A>G\tPASS\t25A>C;25A>G\tPASS\tPanel_v0.2\tV1\n"
            "DPYD\t1\t15\tUNKNOWN\tUNKNOWN\tC\tUNKNOWN\tC\tCAG\trs536\t35A>C;35A>G\tPASS\t35A>C;35A>G?\tUNKNOWN\tPanel_v0.2\tV1\n"
            ".\t15\t24113\tchr15\t684633\tA\tT\tT\tT\trs462;rs9820;rs536\t29482A>T\tPASS\tREF_CALL\tPASS\tPanel_v0.2\tV1\n"
            "BRAF\t2\t154663\tchr2\t40565464\tT\tT\tT\tT\trs154;rs8839\tREF_CALL\tNO_CALL\tREF_CALL\tNO_CALL\tPanel_v0.2\tV1\n"
            "GENE\tX\t15\tchrX\t40\tTT\tAA\tTT\tTT\trs23\tREF_CALL\tNO_CALL\t627AA>TT\tINFERRED_PASS\tPanel_v0.2\tV1\n"
        )
        self.assertEqual(result_expected, result)

    def test_genotype_reporter_non_empty_input_v38(self) -> None:
        all_full_calls = frozenset({
            FullCall(None, None, GeneCoordinate("chr1", 15), "C", ("C", "CAG"),
                     "DPYD", ("rs353",), "3A>C;3A>G?", FullCallFilter.UNKNOWN, "3A>C;3A>G", FullCallFilter.PASS),
            FullCall(GeneCoordinate("1", 5), "A", GeneCoordinate("chr1", 25), "A", ("G", "C"),
                     "DPYD", (".",), "25A>C;25A>G", FullCallFilter.PASS, "25A>C;25A>G", FullCallFilter.PASS),
            FullCall(None, None, GeneCoordinate("chr1", 35), "C", ("C", "CAG"),
                     "DPYD", ("rs536",), "35A>C;35A>G?", FullCallFilter.UNKNOWN, "35A>C;35A>G", FullCallFilter.PASS),
            FullCall(GeneCoordinate("X", 15), "AA", GeneCoordinate("chrX", 40), "TT", ("TT", "TT"),
                     "GENE", ("rs23",), "627AA>TT", FullCallFilter.INFERRED_PASS, "REF_CALL", FullCallFilter.NO_CALL),
            FullCall(GeneCoordinate("2", 154663), "T", GeneCoordinate("chr2", 40565464), "T", ("T", "T"),
                     "BRAF", ("rs154", "rs8839"), "REF_CALL", FullCallFilter.NO_CALL, "REF_CALL", FullCallFilter.NO_CALL),
            FullCall(GeneCoordinate("15", 24113), "T", GeneCoordinate("chr15", 684633), "A", ("T", "T"),
                     ".", ("rs462", "rs9820", "rs536"), "REF_CALL", FullCallFilter.PASS, "29482A>T", FullCallFilter.PASS),
        })
        pgx_analysis = PgxAnalysis(FullCallData(all_full_calls), {})
        panel_id = "Panel_v0.2"
        version = "V1"
        result = GenotypeReporter.get_calls_tsv_text(pgx_analysis, panel_id, version, ReferenceAssembly.V38)

        result_expected = (
            "gene\tchromosome_v37\tposition_v37\tchromosome_v38\tposition_v38\tref_v37\tref_v38\tallele1\tallele2\t"
                "rsid\tvariant_annotation_v37\tfilter_v37\tvariant_annotation_v38\tfilter_v38\tpanel_version\trepo_version\n"
            "DPYD\tUNKNOWN\tUNKNOWN\tchr1\t15\tUNKNOWN\tC\tC\tCAG\trs353\t3A>C;3A>G?\tUNKNOWN\t3A>C;3A>G\tPASS\tPanel_v0.2\tV1\n"
            "DPYD\t1\t5\tchr1\t25\tA\tA\tC\tG\t.\t25A>C;25A>G\tPASS\t25A>C;25A>G\tPASS\tPanel_v0.2\tV1\n"
            "DPYD\tUNKNOWN\tUNKNOWN\tchr1\t35\tUNKNOWN\tC\tC\tCAG\trs536\t35A>C;35A>G?\tUNKNOWN\t35A>C;35A>G\tPASS\tPanel_v0.2\tV1\n"
            ".\t15\t24113\tchr15\t684633\tT\tA\tT\tT\trs462;rs9820;rs536\tREF_CALL\tPASS\t29482A>T\tPASS\tPanel_v0.2\tV1\n"
            "BRAF\t2\t154663\tchr2\t40565464\tT\tT\tT\tT\trs154;rs8839\tREF_CALL\tNO_CALL\tREF_CALL\tNO_CALL\tPanel_v0.2\tV1\n"
            "GENE\tX\t15\tchrX\t40\tAA\tTT\tTT\tTT\trs23\t627AA>TT\tINFERRED_PASS\tREF_CALL\tNO_CALL\tPanel_v0.2\tV1\n"
        )
        self.assertEqual(result_expected, result)

    def test_haplotype_reporter_empty(self) -> None:
        pgx_analysis = PgxAnalysis(FullCallData(frozenset()), {})
        panel = Panel("EmptyPanel", "0.3", frozenset())
        version = "V1"
        result = HaplotypeReporter.get_genotype_tsv_text(pgx_analysis, panel, version)

        result_expected = "gene\thaplotype\tzygosity\tfunction\tlinked_drugs\turl_prescription_info\tpanel_version\trepo_version\n"
        self.assertEqual(result_expected, result)

    def test_haplotype_reporter_non_empty(self) -> None:
        gene_to_haplotype_calls: Dict[str, Set[HaplotypeCall]] = {
            "DPYD": {HaplotypeCall("*2B", 2), HaplotypeCall("*3", 1), HaplotypeCall("*2A", 1)},
            "FAKE": set(),
            "FAKE2": {HaplotypeCall("*1", 1), HaplotypeCall("*4A", 1)},
        }
        pgx_analysis = PgxAnalysis(FullCallData(frozenset()), gene_to_haplotype_calls)
        panel = self.__get_example_panel()
        version = "V1"
        result = HaplotypeReporter.get_genotype_tsv_text(pgx_analysis, panel, version)

        result_expected = (
            "gene\thaplotype\tzygosity\tfunction\tlinked_drugs\turl_prescription_info\tpanel_version\trepo_version\n"
            "DPYD\t*2A\tHET\tNo Function\t5-Fluorouracil;Capecitabine\thttps://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963\tPanel_v0.2\tV1\n"
            "DPYD\t*2B\tHOM\tNo Function\t5-Fluorouracil;Capecitabine\thttps://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963\tPanel_v0.2\tV1\n"
            "DPYD\t*3\tHET\tNormal Function\t5-Fluorouracil;Capecitabine\thttps://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963\tPanel_v0.2\tV1\n"
            "FAKE\tUnresolved Haplotype\tN/A\tUnknown Function\tAspirin\thttps://www.pharmgkb.org/some_other_url\tPanel_v0.2\tV1\n"
            "FAKE2\t*1\tHET\tNormal Function\tAspirin\thttps://www.pharmgkb.org/some_other_url\tPanel_v0.2\tV1\n"
            "FAKE2\t*4A\tHET\tReduced Function\tAspirin\thttps://www.pharmgkb.org/some_other_url\tPanel_v0.2\tV1\n"
        )
        self.assertEqual(result_expected, result)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
