import logging
from typing import NamedTuple, Optional, Tuple

from util.constants import REF_CALL_ANNOTATION_STRING
from util.filter import DualCallFilter
from util.reference_assembly import ReferenceAssembly
from util.reference_site import ReferenceSite
from calls.single_call import AnnotatedSingleCall
from panel.panel import Panel


class Translation(NamedTuple):
    reference_site: Optional[ReferenceSite]
    variant_annotation: Optional[str]
    filter: DualCallFilter


class SingleCallTranslator(object):
    def get_translation(
            self,
            single_call: AnnotatedSingleCall,
            call_reference_assembly: ReferenceAssembly,
            panel: Panel,
    ) -> Translation:
        matching_rs_id = panel.get_perfectly_matching_rs_id(single_call, call_reference_assembly)
        translated_reference_site: Optional[ReferenceSite]
        if matching_rs_id is not None:
            translated_reference_site = panel.get_reference_site(matching_rs_id, call_reference_assembly.opposite())
            has_ref_sequence_difference = panel.has_ref_seq_difference_annotation(matching_rs_id)
            if has_ref_sequence_difference and translated_reference_site is not None:
                result_tuple = self.__get_translated_filter_and_variant_annotation_for_ref_seq_difference(
                    single_call, call_reference_assembly, matching_rs_id, panel, translated_reference_site
                )
                translated_filter, translated_variant_annotation = result_tuple
            elif not has_ref_sequence_difference:
                # no ref seq difference involved
                translated_variant_annotation = single_call.variant_annotation
                if single_call.is_pass():
                    translated_filter = DualCallFilter.PASS
                else:
                    translated_filter = DualCallFilter.NO_CALL
            else:
                error_msg = (
                    f"Cannot handle ref seq difference without knowing the positions vs the other reference assembly: "
                    f"translated_reference_site={translated_reference_site}, call={single_call}"
                )
                raise ValueError(error_msg)
            translation = Translation(translated_reference_site, translated_variant_annotation, translated_filter)
        else:
            # no matching panel variants
            translated_reference_site = None
            if single_call.variant_annotation is None:
                translated_variant_annotation = None
            else:
                translated_variant_annotation = single_call.variant_annotation + "?"
            translated_filter = DualCallFilter.UNKNOWN
            translation = Translation(translated_reference_site, translated_variant_annotation, translated_filter)

            logging.warning(
                f"Unknown variant. Check whether annotation is correct: "
                f"found alleles=({single_call.alleles[0]}, {single_call.alleles[1]}), "
                f"annotation={translated_variant_annotation}"
            )
        return translation

    def __get_translated_filter_and_variant_annotation_for_ref_seq_difference(
            self,
            single_call: AnnotatedSingleCall,
            call_reference_assembly: ReferenceAssembly,
            matching_rs_id: str,
            panel: Panel,
            translated_reference_site: ReferenceSite,
    ) -> Tuple[DualCallFilter, Optional[str]]:
        all_alleles_are_ref_to_opposite_reference_assembly_only = all(
            allele != single_call.reference_site.allele and allele == translated_reference_site.allele
            for allele in single_call.alleles
        )
        all_alleles_are_ref_to_a_reference_assembly = all(
            allele == single_call.reference_site.allele or allele == translated_reference_site.allele
            for allele in single_call.alleles
        )
        translated_variant_annotation: Optional[str]
        if all_alleles_are_ref_to_opposite_reference_assembly_only:
            translated_variant_annotation = REF_CALL_ANNOTATION_STRING
            translated_filter = DualCallFilter.PASS
        elif all_alleles_are_ref_to_a_reference_assembly:
            translated_variant_annotation = panel.get_ref_seq_difference_annotation(
                matching_rs_id, call_reference_assembly.opposite()
            )
            if single_call.is_pass():
                translated_filter = DualCallFilter.PASS
            else:
                translated_filter = DualCallFilter.INFERRED_PASS
        else:
            # unclear what to do
            if single_call.variant_annotation is None:
                translated_variant_annotation = None
            else:
                translated_variant_annotation = single_call.variant_annotation + "?"
            translated_filter = DualCallFilter.UNKNOWN
            logging.warning(
                f"Unexpected allele in ref seq difference location. Check whether annotation is correct: "
                f"found alleles=({single_call.alleles[0]}, {single_call.alleles[1]}), "
                f"annotation={translated_variant_annotation}"
            )
        return translated_filter, translated_variant_annotation
