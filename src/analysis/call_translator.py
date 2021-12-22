import logging
from typing import Set, Tuple, Optional, FrozenSet, NamedTuple

from base.constants import REF_CALL_ANNOTATION_STRING
from base.filter import FullCallFilter, SimpleCallFilter
from base.gene_coordinate import GeneCoordinate
from base.reference_assembly import ReferenceAssembly
from base.reference_site import ReferenceSite
from call_data import SimpleCallData, FullCall, AnnotatedAllele, SimpleCall, FullCallData
from config.panel import Panel


class Translation(NamedTuple):
    reference_site: Optional[ReferenceSite]
    variant_annotation: str
    filter: FullCallFilter


class SimpleCallTranslator(object):
    def get_all_full_call_data(self, simple_call_data: SimpleCallData, panel: Panel) -> FullCallData:
        complete_simple_call_data = self.__add_calls_for_uncalled_variants_in_panel(simple_call_data, panel)
        all_full_calls = self.__get_full_calls_from_simple_calls(complete_simple_call_data, panel)
        return FullCallData(all_full_calls)

    def __add_calls_for_uncalled_variants_in_panel(
            self, simple_call_data: SimpleCallData, panel: Panel
    ) -> SimpleCallData:
        missing_calls = self.__get_calls_for_panel_variants_without_calls(simple_call_data, panel)
        complete_simple_call_data = SimpleCallData(
            simple_call_data.calls.union(missing_calls),
            simple_call_data.reference_assembly,
        )
        return complete_simple_call_data

    def __get_calls_for_panel_variants_without_calls(
            self, simple_call_data: SimpleCallData, panel: Panel
    ) -> FrozenSet[SimpleCall]:
        # assume ref call when no call is found. Set filter to NO_CALL
        reference_assembly = simple_call_data.reference_assembly

        rs_ids_found_in_patient = {rs_id for call in simple_call_data.calls for rs_id in call.rs_ids}
        ref_coordinates_covered_by_found_calls = {
            coordinate
            for call in simple_call_data.calls
            for coordinate in call.reference_site.get_covered_coordinates()
        }

        uncalled_calls = set()
        for gene in panel.get_genes():
            for rs_id in panel.get_rs_ids_for_gene(gene):
                reference_site = panel.get_reference_site(rs_id, reference_assembly)
                relevant_coordinates = reference_site.get_covered_coordinates()
                coordinates_partially_handled = bool(
                    relevant_coordinates.intersection(ref_coordinates_covered_by_found_calls)
                )
                if rs_id not in rs_ids_found_in_patient and not coordinates_partially_handled:
                    # Assuming REF/REF relative to reference assembly

                    uncalled_ref_call = SimpleCall(
                        reference_site,
                        (reference_site.allele, reference_site.allele),
                        gene,
                        (rs_id,),
                        REF_CALL_ANNOTATION_STRING,
                        SimpleCallFilter.NO_CALL,
                    )
                    uncalled_calls.add(uncalled_ref_call)
        return frozenset(uncalled_calls)

    def __get_full_calls_from_simple_calls(self, simple_call_data: SimpleCallData, panel: Panel) -> FrozenSet[FullCall]:
        handled_v37_coordinates: Set[GeneCoordinate] = set()
        handled_v38_coordinates: Set[GeneCoordinate] = set()
        handled_rs_ids: Set[str] = set()
        full_calls = set()
        for simple_call in simple_call_data.calls:
            full_call = self.__get_full_call_from_simple_call(simple_call, panel, simple_call_data.reference_assembly)

            for rs_id in full_call.rs_ids:
                if rs_id in handled_rs_ids:
                    error_msg = (
                        f"Call for rs id that has already been handled:\n"
                        f"call={simple_call}\n"
                        f"handled_rs_ids={handled_rs_ids}"
                    )
                    raise ValueError(error_msg)
                handled_rs_ids.add(rs_id)

            if full_call.reference_site_v37 is not None:
                relevant_v37_coordinates = full_call.reference_site_v37.get_covered_coordinates()
                if relevant_v37_coordinates.intersection(handled_v37_coordinates):
                    warning_msg = (
                        f"Call involves at least one v37 position that has already been handled:\n"
                        f"call={simple_call}\n"
                        f"handled_coords={handled_v37_coordinates}"
                    )
                    logging.warning(warning_msg)
                handled_v37_coordinates.update(relevant_v37_coordinates)
            else:
                warning_msg = f"Could not determine relevant v37 coordinates for call:\ncall={simple_call}"
                logging.warning(warning_msg)

            if full_call.reference_site_v38 is not None:
                relevant_v38_coordinates = full_call.reference_site_v38.get_covered_coordinates()
                if relevant_v38_coordinates.intersection(handled_v38_coordinates):
                    warning_msg = (
                        f"Call involves at least one v38 position that has already been handled:\n"
                        f"call={simple_call}\n"
                        f"handled_coords={handled_v38_coordinates}"
                    )
                    logging.warning(warning_msg)
                handled_v38_coordinates.update(relevant_v38_coordinates)
            else:
                warning_msg = f"Could not determine relevant v38 coordinates for call:\n" f"call={simple_call}"
                logging.warning(warning_msg)

            full_calls.add(full_call)

        return frozenset(full_calls)

    def __get_full_call_from_simple_call(
            self, simple_call: SimpleCall, panel: Panel, call_reference_assembly: ReferenceAssembly
    ) -> FullCall:
        matching_rs_id: Optional[str]
        if panel.contains_rs_id_with_reference_site(simple_call.reference_site, call_reference_assembly):
            matching_rs_id = panel.get_rs_id_with_reference_site(simple_call.reference_site, call_reference_assembly)
        else:
            if any(panel.contains_rs_id(rs_id) for rs_id in simple_call.rs_ids):
                error_msg = f"Rs id is in panel, but the call reference site does not match the panel reference site."
                raise ValueError(error_msg)
            matching_rs_id = None

        complete_simple_call = self.__fill_in_missing_details_from_panel(simple_call, matching_rs_id)
        translation = self.__get_translation(complete_simple_call, call_reference_assembly, matching_rs_id, panel)

        if call_reference_assembly == ReferenceAssembly.V37:
            filter_v37 = self.__get_full_call_filter(complete_simple_call.filter)
            full_call = FullCall(
                reference_site_v37=complete_simple_call.reference_site,
                reference_site_v38=translation.reference_site,
                alleles=complete_simple_call.alleles,
                gene=complete_simple_call.gene,
                rs_ids=complete_simple_call.rs_ids,
                variant_annotation_v37=complete_simple_call.variant_annotation,
                filter_v37=filter_v37,
                variant_annotation_v38=translation.variant_annotation,
                filter_v38=translation.filter,
            )
        elif call_reference_assembly == ReferenceAssembly.V38:
            filter_v38 = self.__get_full_call_filter(complete_simple_call.filter)
            full_call = FullCall(
                reference_site_v37=translation.reference_site,
                reference_site_v38=complete_simple_call.reference_site,
                alleles=complete_simple_call.alleles,
                gene=complete_simple_call.gene,
                rs_ids=complete_simple_call.rs_ids,
                variant_annotation_v37=translation.variant_annotation,
                filter_v37=translation.filter,
                variant_annotation_v38=complete_simple_call.variant_annotation,
                filter_v38=filter_v38,
            )
        else:
            raise NotImplementedError(f"Unrecognized reference assembly: {call_reference_assembly}")

        self.__assert_gene_in_panel(full_call.gene, panel)
        return full_call

    def __get_translation(
            self,
            simple_call: SimpleCall,
            call_reference_assembly: ReferenceAssembly,
            matching_rs_id: Optional[str],
            panel: Panel,
    ) -> Translation:
        translated_reference_site: Optional[ReferenceSite]
        if matching_rs_id is not None:
            translated_reference_site = panel.get_reference_site(matching_rs_id, call_reference_assembly.opposite())
            if panel.has_ref_seq_difference_annotation(matching_rs_id):
                result_tuple = self.__get_translated_filter_and_variant_annotation_for_ref_seq_difference(
                    simple_call, call_reference_assembly, matching_rs_id, panel, translated_reference_site)
                translated_filter, translated_variant_annotation = result_tuple
            else:
                # no ref seq difference involved
                translated_variant_annotation = simple_call.variant_annotation
                if simple_call.is_pass():
                    translated_filter = FullCallFilter.PASS
                else:
                    translated_filter = FullCallFilter.NO_CALL
            translation = Translation(translated_reference_site, translated_variant_annotation, translated_filter)
        else:
            # no matching panel variants
            translated_reference_site = None
            translated_variant_annotation = simple_call.variant_annotation + "?"
            translated_filter = FullCallFilter.UNKNOWN
            translation = Translation(translated_reference_site, translated_variant_annotation, translated_filter)

            logging.warning(
                f"Unknown variant. Check whether annotation is correct: "
                f"found alleles=({simple_call.alleles[0]}, {simple_call.alleles[1]}), "
                f"annotation={translated_variant_annotation}"
            )
        return translation

    def __get_translated_filter_and_variant_annotation_for_ref_seq_difference(
            self,
            simple_call: SimpleCall,
            call_reference_assembly: ReferenceAssembly,
            matching_rs_id: str,
            panel: Panel,
            translated_reference_site: Optional[ReferenceSite],
    ) -> Tuple[FullCallFilter, str]:
        annotated_alleles = self.__get_annotated_alleles(
            simple_call.reference_site,
            simple_call.alleles,
            call_reference_assembly,
            translated_reference_site,
        )
        annotate_as_ref = all(
            self.__is_ref_allele_to_opposite_assembly_due_to_ref_sequence_difference(
                annotated_allele, call_reference_assembly
            )
            for annotated_allele in annotated_alleles
        )
        all_variants_are_ref_to_a_reference_assembly = all(
            self.__allele_is_ref_to_a_reference_assembly(annotated_allele, call_reference_assembly)
            for annotated_allele in annotated_alleles
        )
        if annotate_as_ref:
            translated_variant_annotation = REF_CALL_ANNOTATION_STRING
            translated_filter = FullCallFilter.PASS
        elif all_variants_are_ref_to_a_reference_assembly:
            translated_variant_annotation = panel.get_ref_seq_difference_annotation(
                matching_rs_id, call_reference_assembly.opposite()
            )
            if simple_call.is_pass():
                translated_filter = FullCallFilter.PASS
            else:
                translated_filter = FullCallFilter.INFERRED_PASS
        else:
            translated_variant_annotation = simple_call.variant_annotation + "?"
            translated_filter = FullCallFilter.UNKNOWN
            logging.warning(
                f"Unexpected allele in ref seq difference location. Check whether annotation is correct: "
                f"found alleles=({annotated_alleles[0]}, {annotated_alleles[1]}), "
                f"annotation={translated_variant_annotation}"
            )
        return translated_filter, translated_variant_annotation

    def __fill_in_missing_details_from_panel(
            self, simple_call: SimpleCall, matching_rs_id: Optional[str]
    ) -> SimpleCall:
        if matching_rs_id is not None:
            rs_ids_tuple: Tuple[str, ...]
            if simple_call.rs_ids == tuple():
                rs_ids_tuple = (matching_rs_id,)
            elif simple_call.rs_ids == (matching_rs_id,):
                rs_ids_tuple = simple_call.rs_ids
            else:
                error_msg = (
                    f"Rs id from panel does not match rs ids from VCF: "
                    f"panel_rs_id={matching_rs_id}, vcf_rs_ids={simple_call.rs_ids}"
                )
                raise ValueError(error_msg)

            return SimpleCall(
                simple_call.reference_site,
                simple_call.alleles,
                simple_call.gene,
                rs_ids_tuple,
                simple_call.variant_annotation,
                simple_call.filter,
            )
        else:
            return simple_call

    def __allele_is_ref_to_a_reference_assembly(
            self, annotated_allele: AnnotatedAllele, call_reference_assembly: ReferenceAssembly
    ) -> bool:
        is_ref_vs_call_reference_assembly = not annotated_allele.is_variant_vs(call_reference_assembly)
        is_ref_vs_opposite_reference_assembly = not annotated_allele.is_variant_vs(call_reference_assembly.opposite())
        return is_ref_vs_call_reference_assembly or is_ref_vs_opposite_reference_assembly

    def __is_ref_allele_to_opposite_assembly_due_to_ref_sequence_difference(
            self, annotated_allele: AnnotatedAllele, call_reference_assembly: ReferenceAssembly
    ) -> bool:
        is_variant_vs_call_reference_assembly = annotated_allele.is_variant_vs(call_reference_assembly)
        is_ref_vs_opposite_reference_assembly = not annotated_allele.is_variant_vs(call_reference_assembly.opposite())
        return is_variant_vs_call_reference_assembly and is_ref_vs_opposite_reference_assembly

    def __get_annotated_alleles(
            self,
            call_reference_site: ReferenceSite,
            call_alleles: Tuple[str, str],
            call_reference_assembly: ReferenceAssembly,
            translated_reference_site: Optional[ReferenceSite],
    ) -> Tuple[AnnotatedAllele, AnnotatedAllele]:
        reference_assembly_to_reference_site = {
            call_reference_assembly: call_reference_site,
            call_reference_assembly.opposite(): translated_reference_site,
        }
        annotated_alleles = (
            AnnotatedAllele.from_reference_sites(call_alleles[0], reference_assembly_to_reference_site),
            AnnotatedAllele.from_reference_sites(call_alleles[1], reference_assembly_to_reference_site),
        )
        return annotated_alleles

    def __get_full_call_filter(self, direct_filter: SimpleCallFilter) -> FullCallFilter:
        return FullCallFilter[direct_filter.name]

    def __assert_gene_in_panel(self, gene: str, panel: Panel) -> None:
        if gene not in panel.get_genes():
            error_msg = f"Call for unknown gene:\ngene={gene}"
            raise ValueError(error_msg)
