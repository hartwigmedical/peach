from copy import deepcopy
from typing import Dict, Set, FrozenSet

from call_data import SimpleCallData, FullCall, HaplotypeCall, FullCallData
from call_translator import SimpleCallTranslator
from config.panel import Panel
from haplotype_caller import HaplotypeCaller


class PgxAnalysis(object):
    def __init__(self, full_call_data: FullCallData, gene_to_haplotype_calls: Dict[str, Set[HaplotypeCall]]) -> None:
        self.__full_call_data = full_call_data
        self.__gene_to_haplotype_calls = gene_to_haplotype_calls

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, PgxAnalysis)
            and self.__full_call_data == other.__full_call_data
            and self.__gene_to_haplotype_calls == other.__gene_to_haplotype_calls
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"PgxAnalysis("
            f"full_call_data={self.__full_call_data!r}, "
            f"gene_to_haplotype_calls={self.__gene_to_haplotype_calls!r}, "
            f")"
        )

    def get_all_full_calls(self) -> FrozenSet[FullCall]:
        return self.__full_call_data.calls

    def get_gene_to_haplotype_calls(self) -> Dict[str, Set[HaplotypeCall]]:
        return deepcopy(self.__gene_to_haplotype_calls)


class PgxAnalyser(object):
    @classmethod
    def create_pgx_analysis(cls, vcf_call_data: SimpleCallData, panel: Panel) -> PgxAnalysis:
        full_call_data = SimpleCallTranslator.get_all_full_call_data(vcf_call_data, panel)
        gene_to_haplotype_calls = HaplotypeCaller.get_gene_to_haplotypes_call(full_call_data, panel)
        return PgxAnalysis(full_call_data, gene_to_haplotype_calls)
