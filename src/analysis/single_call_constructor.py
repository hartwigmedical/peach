from typing import Optional, Tuple

from util.reference_assembly import ReferenceAssembly
from calls.single_call import SingleCall
from calls.vcf_call import VcfCall
from panel.panel import Panel


class SingleCallConstructor(object):
    def get_single_call(
            self,
            vcf_call: VcfCall,
            call_reference_assembly: ReferenceAssembly,
            panel: Panel,
    ) -> SingleCall:
        matching_rs_id = panel.get_perfectly_matching_rs_id(vcf_call, call_reference_assembly)
        rs_ids_tuple = self.__get_correct_rs_ids_for_vcf_call(vcf_call, matching_rs_id)
        gene = self.__get_correct_gene_for_vcf_call(vcf_call, call_reference_assembly, panel)

        return SingleCall(
            vcf_call.reference_site,
            vcf_call.alleles,
            gene,
            rs_ids_tuple,
            vcf_call.variant_annotation,
            vcf_call.filter,
        )

    def __get_correct_rs_ids_for_vcf_call(self, vcf_call: VcfCall, matching_rs_id: Optional[str]) -> Tuple[str, ...]:
        rs_ids_tuple: Tuple[str, ...]
        if matching_rs_id is None or vcf_call.rs_ids == (matching_rs_id,):
            rs_ids_tuple = vcf_call.rs_ids
        elif vcf_call.rs_ids == tuple():
            rs_ids_tuple = (matching_rs_id,)
        else:
            error_msg = (
                f"Rs id from panel does not match rs ids from VCF: "
                f"panel_rs_id={matching_rs_id}, vcf_rs_ids={vcf_call.rs_ids}"
            )
            raise ValueError(error_msg)
        return rs_ids_tuple

    def __get_correct_gene_for_vcf_call(
            self, vcf_call: VcfCall, call_reference_assembly: ReferenceAssembly, panel: Panel
    ) -> str:
        relevant_rs_ids_for_call = panel.get_relevant_panel_rs_ids(vcf_call, call_reference_assembly)
        panel_genes_for_call = {panel.get_gene_for_rs_id(rs_id) for rs_id in relevant_rs_ids_for_call}
        if len(panel_genes_for_call) != 1:
            error_msg = (
                f"Call does not match exactly one panel gene: "
                f"panel_genes={sorted(panel_genes_for_call)}, call={vcf_call}"
            )
            raise ValueError(error_msg)
        panel_gene = panel_genes_for_call.pop()
        if vcf_call.gene is not None and panel_gene != vcf_call.gene:
            error_msg = (
                f"Gene name from VCF does not match gene name from panel: "
                f"vcf_gene_name={vcf_call.gene}, panel_gene_name={panel_gene}, call={vcf_call}"
            )
            raise ValueError(error_msg)
        return panel_gene
