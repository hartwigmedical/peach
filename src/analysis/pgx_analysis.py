import logging
from copy import deepcopy
from typing import Dict, Set, FrozenSet

from util.constants import REF_CALL_ANNOTATION_STRING
from util.filter import VcfCallFilter
from calls.haplotype_call import HaplotypeCall
from calls.dual_call import DualCall, DualCallData
from calls.vcf_call import VcfCallData, VcfCall
from analysis.dual_call_constructor import DualCallConstructor
from panel.panel import Panel
from analysis.haplotype_caller import HaplotypeCaller


class PgxAnalysis(object):
    def __init__(self, dual_call_data: DualCallData, gene_to_haplotype_calls: Dict[str, Set[HaplotypeCall]]) -> None:
        self.__dual_call_data = dual_call_data
        self.__gene_to_haplotype_calls = gene_to_haplotype_calls

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, PgxAnalysis)
            and self.__dual_call_data == other.__dual_call_data
            and self.__gene_to_haplotype_calls == other.__gene_to_haplotype_calls
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"PgxAnalysis("
            f"dual_call_data={self.__dual_call_data!r}, "
            f"gene_to_haplotype_calls={self.__gene_to_haplotype_calls!r}, "
            f")"
        )

    def get_all_dual_calls(self) -> FrozenSet[DualCall]:
        return self.__dual_call_data.calls

    def get_gene_to_haplotype_calls(self) -> Dict[str, Set[HaplotypeCall]]:
        return deepcopy(self.__gene_to_haplotype_calls)


class PgxAnalyser(object):
    def create_pgx_analysis(self, vcf_call_data: VcfCallData, panel: Panel) -> PgxAnalysis:
        logging.info("Adding calls for panel variants without calls")
        complete_vcf_call_data = self.__add_calls_for_uncalled_variants_in_panel(vcf_call_data, panel)

        logging.info(f"Annotating call data for reference assembly {vcf_call_data.reference_assembly.opposite().name}")
        dual_call_data = DualCallConstructor().get_dual_call_data(complete_vcf_call_data, panel)

        logging.info(f"Calling haplotypes")
        gene_to_haplotype_calls = HaplotypeCaller().get_gene_to_haplotypes_call(dual_call_data, panel)

        return PgxAnalysis(dual_call_data, gene_to_haplotype_calls)

    def __add_calls_for_uncalled_variants_in_panel(
            self, vcf_call_data: VcfCallData, panel: Panel
    ) -> VcfCallData:
        missing_calls = self.__get_calls_for_panel_variants_without_calls(vcf_call_data, panel)
        complete_vcf_call_data = VcfCallData(
            vcf_call_data.calls.union(missing_calls),
            vcf_call_data.reference_assembly,
        )
        return complete_vcf_call_data

    def __get_calls_for_panel_variants_without_calls(
            self, vcf_call_data: VcfCallData, panel: Panel
    ) -> FrozenSet[VcfCall]:
        # assume ref call when no call is found. Set filter to NO_CALL
        rs_ids_found_in_patient = {rs_id for call in vcf_call_data.calls for rs_id in call.rs_ids}
        ref_coordinates_covered_by_found_calls = {
            coordinate
            for call in vcf_call_data.calls
            for coordinate in call.reference_site.get_covered_coordinates()
        }

        uncalled_calls = set()
        for gene in panel.get_genes():
            for rs_id in panel.get_rs_ids_for_gene(gene):
                reference_site = panel.get_reference_site(rs_id, vcf_call_data.reference_assembly)
                relevant_coordinates = reference_site.get_covered_coordinates()
                coordinates_partially_handled = bool(
                    relevant_coordinates.intersection(ref_coordinates_covered_by_found_calls)
                )
                if rs_id not in rs_ids_found_in_patient and not coordinates_partially_handled:
                    # Assuming REF/REF relative to reference assembly

                    uncalled_ref_call = VcfCall(
                        reference_site,
                        (reference_site.allele, reference_site.allele),
                        gene,
                        (rs_id,),
                        REF_CALL_ANNOTATION_STRING,
                        VcfCallFilter.NO_CALL,
                    )
                    uncalled_calls.add(uncalled_ref_call)
        return frozenset(uncalled_calls)
