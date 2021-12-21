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
    @classmethod
    def get_all_full_call_data(cls, simple_call_data: SimpleCallData, panel: Panel) -> FullCallData:
        complete_simple_call_data = cls.__add_calls_for_uncalled_variants_in_panel(simple_call_data, panel)
        all_full_calls = cls.__get_full_calls_from_simple_calls(complete_simple_call_data, panel)
        return FullCallData(all_full_calls)

    @classmethod
    def __add_calls_for_uncalled_variants_in_panel(
        cls, simple_call_data: SimpleCallData, panel: Panel
    ) -> SimpleCallData:
        missing_calls = cls.__get_calls_for_panel_variants_without_calls(simple_call_data, panel)
        complete_simple_call_data = SimpleCallData(
            simple_call_data.calls.union(missing_calls),
            simple_call_data.reference_assembly,
        )
        return complete_simple_call_data

    @classmethod
    def __get_calls_for_panel_variants_without_calls(
        cls, simple_call_data: SimpleCallData, panel: Panel
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
        for gene_info in panel.get_gene_infos():
            for rs_id_info in gene_info.rs_id_infos:
                relevant_coordinates = rs_id_info.get_reference_site(reference_assembly).get_covered_coordinates()
                coordinates_partially_handled = bool(
                    relevant_coordinates.intersection(ref_coordinates_covered_by_found_calls)
                )
                if rs_id_info.rs_id not in rs_ids_found_in_patient and not coordinates_partially_handled:
                    # Assuming REF/REF relative to reference assembly
                    reference_site = rs_id_info.get_reference_site(reference_assembly)
                    uncalled_ref_call = SimpleCall(
                        reference_site,
                        (reference_site.allele, reference_site.allele),
                        gene_info.gene,
                        (rs_id_info.rs_id,),
                        REF_CALL_ANNOTATION_STRING,
                        SimpleCallFilter.NO_CALL,
                    )
                    uncalled_calls.add(uncalled_ref_call)
        return frozenset(uncalled_calls)

    @classmethod
    def __get_full_calls_from_simple_calls(cls, simple_call_data: SimpleCallData, panel: Panel) -> FrozenSet[FullCall]:
        handled_v37_coordinates: Set[GeneCoordinate] = set()
        handled_v38_coordinates: Set[GeneCoordinate] = set()
        handled_rs_ids: Set[str] = set()
        full_calls = set()
        for simple_call in simple_call_data.calls:
            full_call = cls.__get_full_call_from_simple_call(simple_call, panel, simple_call_data.reference_assembly)

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

    @classmethod
    def __get_full_call_from_simple_call(
        cls, simple_call: SimpleCall, panel: Panel, call_reference_assembly: ReferenceAssembly
    ) -> FullCall:
        cls.__assert_gene_in_panel(simple_call.gene, panel)

        complete_simple_call = cls.__fill_in_rs_ids_if_needed(simple_call, panel, call_reference_assembly)
        translation = cls.__get_translation_to_other_assembly(complete_simple_call, panel, call_reference_assembly)

        if call_reference_assembly == ReferenceAssembly.V37:
            filter_v37 = cls.__get_full_call_filter(complete_simple_call.filter)
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
            filter_v38 = cls.__get_full_call_filter(complete_simple_call.filter)
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
        return full_call

    @classmethod
    def __fill_in_rs_ids_if_needed(
        cls, call: SimpleCall, panel: Panel, reference_assembly: ReferenceAssembly
    ) -> SimpleCall:
        rs_ids: Tuple[str, ...]
        if call.rs_ids == tuple() and panel.contains_rs_id_matching_call(call, reference_assembly):
            rs_id_info = panel.get_matching_rs_id_info(call.reference_site, reference_assembly)
            rs_ids = (rs_id_info.rs_id,)
            new_simple_call = SimpleCall(
                call.reference_site,
                call.alleles,
                call.gene,
                rs_ids,
                call.variant_annotation,
                call.filter,
            )
            return new_simple_call
        else:
            return call

    @classmethod
    def __get_translation_to_other_assembly(
        cls, call: SimpleCall, panel: Panel, call_reference_assembly: ReferenceAssembly
    ) -> Translation:
        translated_reference_site = cls.__get_translated_reference_site(call, panel, call_reference_assembly)
        translated_filter, translated_variant_annotation = cls.__get_translated_filter_and_variant_annotation(
            call, panel, call_reference_assembly, translated_reference_site
        )
        translation = Translation(translated_reference_site, translated_variant_annotation, translated_filter)
        return translation

    @classmethod
    def __get_translated_reference_site(
        cls, call: SimpleCall, panel: Panel, call_reference_assembly: ReferenceAssembly
    ) -> Optional[ReferenceSite]:
        if panel.contains_rs_id_matching_call(call, call_reference_assembly):
            rs_id_info = panel.get_matching_rs_id_info(call.reference_site, call_reference_assembly)
            cls.__assert_rs_id_call_matches_info(call.rs_ids, (rs_id_info.rs_id,))
            return rs_id_info.get_reference_site(call_reference_assembly.opposite())
        else:
            # unknown variant
            return None

    @classmethod
    def __get_translated_filter_and_variant_annotation(
        cls,
        call: SimpleCall,
        panel: Panel,
        call_reference_assembly: ReferenceAssembly,
        translated_reference_site: Optional[ReferenceSite],
    ) -> Tuple[FullCallFilter, str]:
        annotated_alleles = cls.__get_annotated_alleles(call, call_reference_assembly, translated_reference_site)
        if panel.has_ref_seq_difference_annotation(call.gene, call.reference_site, call_reference_assembly):
            annotate_as_ref = all(
                cls.__is_ref_allele_to_opposite_assembly_due_to_ref_sequence_difference(
                    annotated_allele, call_reference_assembly
                )
                for annotated_allele in annotated_alleles
            )
            all_variants_are_ref_to_a_reference_assembly = all(
                cls.__allele_is_ref_to_a_reference_assembly(annotated_allele, call_reference_assembly)
                for annotated_allele in annotated_alleles
            )
            if annotate_as_ref:
                translated_variant_annotation = REF_CALL_ANNOTATION_STRING
                translated_filter = FullCallFilter.PASS
            elif all_variants_are_ref_to_a_reference_assembly:
                translated_variant_annotation = panel.get_ref_seq_difference_annotation(
                    call.gene, call.reference_site, call_reference_assembly
                )
                if call.is_pass():
                    translated_filter = FullCallFilter.PASS
                else:
                    translated_filter = FullCallFilter.INFERRED_PASS
            else:
                translated_variant_annotation = call.variant_annotation + "?"
                translated_filter = FullCallFilter.UNKNOWN
                logging.warning(
                    f"Unexpected allele in ref seq difference location. Check whether annotation is correct: "
                    f"found alleles=({annotated_alleles[0]}, {annotated_alleles[1]}), "
                    f"annotation={translated_variant_annotation}"
                )
        elif panel.contains_rs_id_matching_call(call, call_reference_assembly):
            # known variant and no ref seq differences involved
            translated_variant_annotation = call.variant_annotation
            if call.is_pass():
                translated_filter = FullCallFilter.PASS
            else:
                translated_filter = FullCallFilter.NO_CALL
        else:
            # unknown variant, no ref seq difference involved
            translated_variant_annotation = call.variant_annotation + "?"
            translated_filter = FullCallFilter.UNKNOWN
            logging.warning(
                f"Unknown variant. Check whether annotation is correct: "
                f"found alleles=({annotated_alleles[0]}, {annotated_alleles[1]}), "
                f"annotation={translated_variant_annotation}"
            )
        return translated_filter, translated_variant_annotation

    @classmethod
    def __allele_is_ref_to_a_reference_assembly(
        cls, annotated_allele: AnnotatedAllele, call_reference_assembly: ReferenceAssembly
    ) -> bool:
        is_ref_vs_call_reference_assembly = not annotated_allele.is_variant_vs(call_reference_assembly)
        is_ref_vs_opposite_reference_assembly = not annotated_allele.is_variant_vs(call_reference_assembly.opposite())
        return is_ref_vs_call_reference_assembly or is_ref_vs_opposite_reference_assembly

    @classmethod
    def __is_ref_allele_to_opposite_assembly_due_to_ref_sequence_difference(
        cls, annotated_allele: AnnotatedAllele, call_reference_assembly: ReferenceAssembly
    ) -> bool:
        is_variant_vs_call_reference_assembly = annotated_allele.is_variant_vs(call_reference_assembly)
        is_ref_vs_opposite_reference_assembly = not annotated_allele.is_variant_vs(call_reference_assembly.opposite())
        return is_variant_vs_call_reference_assembly and is_ref_vs_opposite_reference_assembly

    @classmethod
    def __get_annotated_alleles(
        cls,
        call: SimpleCall,
        call_reference_assembly: ReferenceAssembly,
        translated_reference_site: Optional[ReferenceSite],
    ) -> Tuple[AnnotatedAllele, AnnotatedAllele]:
        reference_assembly_to_reference_site = {
            call_reference_assembly: call.reference_site,
            call_reference_assembly.opposite(): translated_reference_site,
        }
        annotated_alleles = (
            AnnotatedAllele.from_reference_sites(call.alleles[0], reference_assembly_to_reference_site),
            AnnotatedAllele.from_reference_sites(call.alleles[1], reference_assembly_to_reference_site),
        )
        return annotated_alleles

    @classmethod
    def __get_full_call_filter(cls, direct_filter: SimpleCallFilter) -> FullCallFilter:
        return FullCallFilter[direct_filter.name]

    @classmethod
    def __assert_rs_id_call_matches_info(cls, rs_ids_call: Tuple[str, ...], rs_ids_info: Tuple[str, ...]) -> None:
        if rs_ids_call != tuple() and rs_ids_call != rs_ids_info:
            # TODO: make this more flexible, if necessary
            error_msg = (
                f"Given rs id does not match rs id from panel: from call={rs_ids_call}, from panel={rs_ids_info}"
            )
            raise ValueError(error_msg)

    @classmethod
    def __assert_gene_in_panel(cls, gene: str, panel: Panel) -> None:
        if gene not in panel.get_genes():
            error_msg = f"Call for unknown gene:\ngene={gene}"
            raise ValueError(error_msg)
