from typing import Set, Tuple, Optional, FrozenSet, NamedTuple

from base.constants import REF_CALL_ANNOTATION_STRING
from base.filter import FullCallFilter, SimpleCallFilter
from base.gene_coordinate import GeneCoordinate
from base.reference_assembly import ReferenceAssembly
from call_data import SimpleCallData, FullCall, AnnotatedAllele, SimpleCall, FullCallData
from config.panel import Panel


class Translation(NamedTuple):
    start_coordinate: Optional[GeneCoordinate]
    reference_allele: Optional[str]
    variant_annotation: str
    filter: FullCallFilter


class V37CallTranslator(object):
    @classmethod
    def get_all_full_call_data(cls, simple_call_data: SimpleCallData, panel: Panel) -> FullCallData:
        complete_simple_call_data = cls.__add_calls_for_uncalled_variants_in_panel(simple_call_data, panel)
        if simple_call_data.reference_assembly != ReferenceAssembly.V37:
            raise NotImplementedError("WIP")
        all_full_calls = cls.__get_full_calls_from_simple_calls(complete_simple_call_data, panel)
        return FullCallData(all_full_calls)

    @classmethod
    def __add_calls_for_uncalled_variants_in_panel(cls, simple_call_data: SimpleCallData, panel: Panel) -> SimpleCallData:
        missing_calls = cls.__get_calls_for_panel_variants_without_calls(simple_call_data, panel)
        complete_simple_call_data = SimpleCallData(
            simple_call_data.calls.union(missing_calls),
            simple_call_data.reference_assembly,
        )
        return complete_simple_call_data

    @classmethod
    def __get_calls_for_panel_variants_without_calls(
            cls, simple_call_data: SimpleCallData, panel: Panel) -> FrozenSet[SimpleCall]:
        # assume ref call when no call is found. Set filter to NO_CALL
        reference_assembly = simple_call_data.reference_assembly
        
        rs_ids_found_in_patient = {
            rs_id for call in simple_call_data.calls for rs_id in call.rs_ids if rs_id != "."
        }
        coordinates_covered_by_found_calls = {
            coordinate for call in simple_call_data.calls for coordinate in call.get_relevant_coordinates()
        }

        uncalled_calls = set()
        for gene_info in panel.get_gene_infos():
            for rs_id_info in gene_info.rs_id_infos:
                coordinates_partially_handled = bool(
                    rs_id_info.get_relevant_coordinates(reference_assembly).intersection(
                        coordinates_covered_by_found_calls)
                )
                if rs_id_info.rs_id not in rs_ids_found_in_patient and not coordinates_partially_handled:
                    # Assuming REF/REF relative to reference assembly
                    start_coordinate = rs_id_info.get_start_coordinate(reference_assembly)
                    reference_allele = rs_id_info.get_reference_allele(reference_assembly)
                    uncalled_ref_call = SimpleCall(
                        start_coordinate,
                        reference_allele,
                        (reference_allele, reference_allele),
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
            full_call = cls.__get_full_call_from_v37_call(simple_call, panel, simple_call_data.reference_assembly)

            for rs_id in full_call.rs_ids:
                if rs_id != ".":
                    if rs_id in handled_rs_ids:
                        error_msg = (
                            f"Call for rs id that has already been handled:\n"
                            f"call={simple_call}\n"
                            f"handled_rs_ids={handled_rs_ids}")
                        raise ValueError(error_msg)
                    handled_rs_ids.add(rs_id)

            relevant_v37_coordinates = full_call.get_relevant_v37_coordinates()
            if relevant_v37_coordinates is not None:
                if relevant_v37_coordinates.intersection(handled_v37_coordinates):
                    warning_msg = (
                        f"[WARN] Call involves at least one v37 position that has already been handled:\n"
                        f"call={simple_call}\n"
                        f"handled_coords={handled_v37_coordinates}")
                    print(warning_msg)
                handled_v37_coordinates.update(relevant_v37_coordinates)
            else:
                warning_msg = (
                    f"[WARN] Could not determine relevant v37 coordinates for call:\n"
                    f"call={simple_call}")
                print(warning_msg)

            relevant_v38_coordinates = full_call.get_relevant_v38_coordinates()
            if relevant_v38_coordinates is not None:
                if relevant_v38_coordinates.intersection(handled_v38_coordinates):
                    warning_msg = (
                        f"[WARN] Call involves at least one v38 position that has already been handled:\n"
                        f"call={simple_call}\n"
                        f"handled_coords={handled_v38_coordinates}")
                    print(warning_msg)
                handled_v38_coordinates.update(relevant_v38_coordinates)
            else:
                warning_msg = (
                    f"[WARN] Could not determine relevant v38 coordinates for call:\n"
                    f"call={simple_call}")
                print(warning_msg)

            full_calls.add(full_call)

        return frozenset(full_calls)

    @classmethod
    def __get_full_call_from_v37_call(
            cls, v37_call: SimpleCall, panel: Panel, call_reference_assembly: ReferenceAssembly) -> FullCall:
        cls.__assert_gene_in_panel(v37_call.gene, panel)

        v37_call = cls.__fill_in_rs_ids_if_needed(v37_call, panel, call_reference_assembly)
        translation = cls.get_translation_to_other_assembly(v37_call, panel, call_reference_assembly)
        filter_v37 = cls.__get_full_call_filter(v37_call.filter)
        full_call = FullCall(
            v37_call.start_coordinate, v37_call.reference_allele,
            translation.start_coordinate, translation.reference_allele,
            v37_call.alleles, v37_call.gene, v37_call.rs_ids,
            v37_call.variant_annotation, filter_v37,
            translation.variant_annotation, translation.filter,
        )
        return full_call

    @classmethod
    def __fill_in_rs_ids_if_needed(
            cls, call: SimpleCall, panel: Panel, reference_assembly: ReferenceAssembly) -> SimpleCall:
        rs_ids: Tuple[str, ...]
        if call.rs_ids == (".",) and panel.contains_rs_id_matching_call(call, reference_assembly):
            rs_id_info = panel.get_matching_rs_id_info(call.start_coordinate, call.reference_allele, reference_assembly)
            rs_ids = (rs_id_info.rs_id,)
            new_simple_call = SimpleCall(
                call.start_coordinate,
                call.reference_allele,
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
    def get_translation_to_other_assembly(
            cls, call: SimpleCall, panel: Panel, call_reference_assembly: ReferenceAssembly) -> Translation:
        # determine start_coordinate_v38, reference_allele_v38
        translated_start_coordinate: Optional[GeneCoordinate]
        translated_reference_allele: Optional[str]
        if panel.contains_rs_id_matching_call(call, call_reference_assembly):
            rs_id_info = panel.get_matching_rs_id_info(
                call.start_coordinate, call.reference_allele, call_reference_assembly)
            cls.__assert_rs_id_call_matches_info(call.rs_ids, (rs_id_info.rs_id,))

            translated_start_coordinate = rs_id_info.get_start_coordinate(call_reference_assembly.opposite())
            translated_reference_allele = rs_id_info.get_reference_allele(call_reference_assembly.opposite())
        else:
            # unknown variant
            translated_start_coordinate = None
            translated_reference_allele = None

        if call_reference_assembly != ReferenceAssembly.V37:
            raise NotImplementedError("WIP")
        
        annotated_alleles = (
            AnnotatedAllele.from_alleles(call.alleles[0], call.reference_allele, translated_reference_allele),
            AnnotatedAllele.from_alleles(call.alleles[1], call.reference_allele, translated_reference_allele),
        )
        # determine variant annotation v38 and filter v38
        if panel.has_ref_seq_difference_annotation(
                call.gene, call.start_coordinate, call.reference_allele, call_reference_assembly):
            v38_ref_call_due_to_ref_sequence_difference = all(
                annotated.is_annotated_vs_v37()
                and annotated.is_annotated_vs_v38()
                and annotated.is_variant_vs_v37
                and not annotated.is_variant_vs_v38
                for annotated in annotated_alleles
            )
            all_variants_ref_to_v37_or_v38 = all(
                (annotated.is_annotated_vs_v37() and not annotated.is_variant_vs_v37)
                or (annotated.is_annotated_vs_v38() and not annotated.is_variant_vs_v38)
                for annotated in annotated_alleles
            )
            if v38_ref_call_due_to_ref_sequence_difference:
                variant_annotation_v38 = REF_CALL_ANNOTATION_STRING
                filter_v38 = FullCallFilter.PASS
            elif all_variants_ref_to_v37_or_v38:
                variant_annotation_v38 = panel.get_ref_seq_difference_annotation(
                    call.gene, call.start_coordinate, call.reference_allele, call_reference_assembly)
                if call.is_pass():
                    filter_v38 = FullCallFilter.PASS
                else:
                    filter_v38 = FullCallFilter.INFERRED_PASS
            else:
                variant_annotation_v38 = call.variant_annotation + "?"
                filter_v38 = FullCallFilter.UNKNOWN
                print(
                    f"[WARN] Unexpected allele in ref seq difference location. Check whether annotation is correct: "
                    f"found alleles=({annotated_alleles[0]}, {annotated_alleles[1]}), "
                    f"annotation={variant_annotation_v38}"
                )
        elif panel.contains_rs_id_matching_call(call, call_reference_assembly):
            # known variant and no ref seq differences involved
            variant_annotation_v38 = call.variant_annotation
            if call.is_pass():
                filter_v38 = FullCallFilter.PASS
            else:
                filter_v38 = FullCallFilter.NO_CALL
        else:
            # unknown variant, no ref seq difference involved
            variant_annotation_v38 = call.variant_annotation + "?"
            filter_v38 = FullCallFilter.UNKNOWN
            print(
                f"[WARN] Unknown variant. Check whether annotation is correct: "
                f"found alleles=({annotated_alleles[0]}, {annotated_alleles[1]}), "
                f"annotation={variant_annotation_v38}"
            )
        translation = Translation(translated_start_coordinate, translated_reference_allele, variant_annotation_v38, filter_v38)
        return translation

    @classmethod
    def __get_full_call_filter(cls, direct_filter: SimpleCallFilter) -> FullCallFilter:
        return FullCallFilter[direct_filter.name]

    @classmethod
    def __assert_rs_id_call_matches_info(cls, rs_ids_call: Tuple[str, ...], rs_ids_info: Tuple[str, ...]) -> None:
        if rs_ids_call != (".",) and rs_ids_call != rs_ids_info:
            # TODO: make this more flexible, if necessary
            error_msg = (f"Given rs id does not match rs id from panel: "
                         f"from call={rs_ids_call}, from panel={rs_ids_info}")
            raise ValueError(error_msg)

    @classmethod
    def __assert_gene_in_panel(cls, gene: str, panel: Panel) -> None:
        if gene not in panel.get_genes():
            error_msg = f"Call for unknown gene:\ngene={gene}"
            raise ValueError(error_msg)
