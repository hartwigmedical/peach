import collections
import logging
from copy import deepcopy
from typing import Dict, Set, DefaultDict, Tuple

from util.gene_coordinate import GeneCoordinate
from util.reference_assembly import ReferenceAssembly
from calls.haplotype_call import HaplotypeCall
from calls.dual_call import DualCall, DualCallData
from panel.panel import Panel
from panel.variant import Variant


class HaplotypeCaller(object):
    HAPLOTYPE_CALLING_REFERENCE_ASSEMBLY = ReferenceAssembly.V38

    def get_gene_to_haplotypes_call(self, dual_call_data: DualCallData, panel: Panel) -> Dict[str, Set[HaplotypeCall]]:
        self.__assert_call_data_is_handleable(dual_call_data)

        gene_to_haplotype_calls = {}
        for gene in panel.get_genes():
            logging.info(f"Calling haplotypes for {gene}")
            gene_to_haplotype_calls[gene] = self.__get_haplotypes_call(gene, dual_call_data, panel)
        return gene_to_haplotype_calls

    def __get_haplotypes_call(self, gene: str, dual_call_data: DualCallData, panel: Panel) -> Set[HaplotypeCall]:
        try:
            variant_to_count = self.__get_variant_to_count_for_gene(dual_call_data, gene)

            explaining_haplotype_combinations = self.__get_explaining_haplotype_combinations(
                variant_to_count, gene, panel
            )

            if not explaining_haplotype_combinations:
                error_msg = f"No explaining haplotype combinations"
                raise ValueError(error_msg)

            minimal_explaining_haplotype_combination = self.__get_minimal_haplotype_combination(
                explaining_haplotype_combinations
            )

            haplotype_calls = self.__get_haplotype_calls_from_haplotype_names(
                minimal_explaining_haplotype_combination, panel.get_wild_type_haplotype_name(gene)
            )

            return haplotype_calls

        except Exception as e:
            logging.error(f"Cannot resolve haplotype for gene {gene}. Error: {e}")
            return set()

    def __get_variant_to_count_for_gene(self, dual_call_data: DualCallData, gene: str) -> DefaultDict[Variant, int]:
        dual_calls_for_gene = {call for call in dual_call_data.calls if call.gene == gene}
        variant_to_count: DefaultDict[Variant, int] = collections.defaultdict(int)
        for call in dual_calls_for_gene:
            self.__assert_handleable_call(call)
            reference_site = call.get_reference_site(self.HAPLOTYPE_CALLING_REFERENCE_ASSEMBLY)
            if reference_site is None:
                error_msg = (
                    f"Cannot call haplotypes, since call is "
                    f"not annotated vs {self.HAPLOTYPE_CALLING_REFERENCE_ASSEMBLY}: call={call}"
                )
                raise ValueError(error_msg)
            for allele in call.alleles:
                if allele != reference_site.allele:
                    variant_to_count[Variant(call.rs_ids[0], allele)] += 1
        return variant_to_count

    def __get_explaining_haplotype_combinations(
        self, variant_to_count: DefaultDict[Variant, int], gene: str, panel: Panel
    ) -> Set[Tuple[str, ...]]:
        """
        Gets combinations of haplotypes that explain all variants in the stated amounts. Uses recursion.
        Always makes sure that the haplotypes in a haplotype combination are ordered alphabetically to
        ensure that each haplotype combination exists only once in the result set.
        """
        if any(count < 0 for count in variant_to_count.values()):
            return set()
        if all(count == 0 for count in variant_to_count.values()):
            return {tuple()}

        result_set = set()
        for haplotype_name in panel.get_non_wild_type_haplotype_names(gene):
            reduced_variant_to_count = deepcopy(variant_to_count)
            for variant in panel.get_variants_for_haplotype(gene, haplotype_name):
                reduced_variant_to_count[variant] -= 1

            combinations_for_reduced_variant_set = self.__get_explaining_haplotype_combinations(
                reduced_variant_to_count, gene, panel
            )
            for combination in combinations_for_reduced_variant_set:
                result_set.add(tuple(sorted(list(combination) + [haplotype_name])))

        return result_set

    def __get_minimal_haplotype_combination(
        self, explaining_haplotype_combinations: Set[Tuple[str, ...]]
    ) -> Tuple[str, ...]:
        min_haplotype_count = min(len(combination) for combination in explaining_haplotype_combinations)
        minimal_explaining_haplotype_combinations = {
            combination for combination in explaining_haplotype_combinations if len(combination) == min_haplotype_count
        }
        if len(minimal_explaining_haplotype_combinations) > 1:
            error_msg = (
                f"No unique minimal explaining haplotype combination: "
                f"options={minimal_explaining_haplotype_combinations}"
            )
            raise ValueError(error_msg)
        minimal_explaining_haplotype_combination = minimal_explaining_haplotype_combinations.pop()
        return minimal_explaining_haplotype_combination

    def __get_haplotype_calls_from_haplotype_names(
        self, haplotype_name_combination: Tuple[str, ...], wild_type_haplotype_name: str
    ) -> Set[HaplotypeCall]:
        haplotype_to_count: DefaultDict[str, int] = collections.defaultdict(int)
        for haplotype in haplotype_name_combination:
            haplotype_to_count[haplotype] += 1

        haplotype_calls = set()
        for haplotype, count in haplotype_to_count.items():
            if count == 1 or count == 2:
                haplotype_calls.add(HaplotypeCall(haplotype, count))
            else:
                error_msg = f"Impossible count for haplotype: haplotype={haplotype}, count={count}"
                raise ValueError(error_msg)

        called_haplotypes_count = sum(haplotype_to_count.values())
        if called_haplotypes_count == 0:
            haplotype_calls.add(HaplotypeCall(wild_type_haplotype_name, 2))
        elif called_haplotypes_count == 1:
            haplotype_calls.add(HaplotypeCall(wild_type_haplotype_name, 1))

        return haplotype_calls

    def __assert_handleable_call(self, call: DualCall) -> None:
        if len(call.rs_ids) > 1:
            error_msg = f"Call has more than one rs id: rs ids={call.rs_ids}, call={call}"
            raise ValueError(error_msg)
        if len(call.rs_ids) < 1:
            error_msg = f"Call has zero rs ids: call={call}"
            raise ValueError(error_msg)
        if call.rs_ids == tuple():
            error_msg = f"Call has unknown rs id: call={call}"
            raise ValueError(error_msg)

    def __assert_call_data_is_handleable(self, dual_call_data: DualCallData) -> None:
        handled_v37_coordinates: Set[GeneCoordinate] = set()
        handled_v38_coordinates: Set[GeneCoordinate] = set()
        handled_rs_ids: Set[str] = set()
        for dual_call in dual_call_data.calls:
            for rs_id in dual_call.rs_ids:
                if rs_id in handled_rs_ids:
                    error_msg = (
                        f"Call for rs id that has already been handled:\n"
                        f"call={dual_call}\n"
                        f"handled_rs_ids={handled_rs_ids}"
                    )
                    raise ValueError(error_msg)
                handled_rs_ids.add(rs_id)

            if dual_call.reference_site_v37 is not None:
                relevant_v37_coordinates = dual_call.reference_site_v37.get_covered_coordinates()
                if relevant_v37_coordinates.intersection(handled_v37_coordinates):
                    warning_msg = (
                        f"Call involves at least one v37 position that has already been handled:\n"
                        f"call={dual_call}\n"
                        f"handled_coordinates={handled_v37_coordinates}"
                    )
                    logging.warning(warning_msg)
                handled_v37_coordinates.update(relevant_v37_coordinates)
            else:
                warning_msg = f"Could not determine relevant v37 coordinates for call:\ncall={dual_call}"
                logging.warning(warning_msg)

            if dual_call.reference_site_v38 is not None:
                relevant_v38_coordinates = dual_call.reference_site_v38.get_covered_coordinates()
                if relevant_v38_coordinates.intersection(handled_v38_coordinates):
                    warning_msg = (
                        f"Call involves at least one v38 position that has already been handled:\n"
                        f"call={dual_call}\n"
                        f"handled_coordinates={handled_v38_coordinates}"
                    )
                    logging.warning(warning_msg)
                handled_v38_coordinates.update(relevant_v38_coordinates)
            else:
                warning_msg = f"Could not determine relevant v38 coordinates for call:\n" f"call={dual_call}"
                logging.warning(warning_msg)

