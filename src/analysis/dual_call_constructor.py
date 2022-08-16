from analysis.single_call_constructor import SingleCallConstructor
from analysis.single_call_translator import SingleCallTranslator
from util.filter import DualCallFilter, VcfCallFilter
from util.reference_assembly import ReferenceAssembly
from calls.dual_call import DualCall, DualCallData
from calls.vcf_call import VcfCall, VcfCallData
from panel.panel import Panel


class DualCallConstructor(object):
    def get_dual_call_data(self, vcf_call_data: VcfCallData, panel: Panel) -> DualCallData:
        dual_calls = set()
        for vcf_call in vcf_call_data.calls:
            dual_call = self.__get_dual_call_from_vcf_call(vcf_call, vcf_call_data.reference_assembly, panel)
            dual_calls.add(dual_call)

        return DualCallData(frozenset(dual_calls))

    def __get_dual_call_from_vcf_call(
            self, vcf_call: VcfCall, call_reference_assembly: ReferenceAssembly, panel: Panel
    ) -> DualCall:
        single_call = SingleCallConstructor().get_single_call(vcf_call, call_reference_assembly, panel)
        translation = SingleCallTranslator().get_translation(single_call, call_reference_assembly, panel)

        if call_reference_assembly == ReferenceAssembly.V37:
            filter_v37 = self.__get_dual_call_filter(single_call.filter)
            dual_call = DualCall(
                reference_site_v37=single_call.reference_site,
                reference_site_v38=translation.reference_site,
                alleles=single_call.alleles,
                gene=single_call.gene,
                rs_ids=single_call.rs_ids,
                variant_annotation_v37=single_call.variant_annotation,
                filter_v37=filter_v37,
                variant_annotation_v38=translation.variant_annotation,
                filter_v38=translation.filter,
            )
        elif call_reference_assembly == ReferenceAssembly.V38:
            filter_v38 = self.__get_dual_call_filter(single_call.filter)
            dual_call = DualCall(
                reference_site_v37=translation.reference_site,
                reference_site_v38=single_call.reference_site,
                alleles=single_call.alleles,
                gene=single_call.gene,
                rs_ids=single_call.rs_ids,
                variant_annotation_v37=translation.variant_annotation,
                filter_v37=translation.filter,
                variant_annotation_v38=single_call.variant_annotation,
                filter_v38=filter_v38,
            )
        else:
            raise NotImplementedError(f"Unrecognized reference assembly: {call_reference_assembly}")

        return dual_call

    def __get_dual_call_filter(self, direct_filter: VcfCallFilter) -> DualCallFilter:
        return DualCallFilter[direct_filter.name]
