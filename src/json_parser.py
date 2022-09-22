from typing import Optional, Any, Dict, FrozenSet

from calls.dual_call import DualCall
from util.gene_coordinate import GeneCoordinate
from util.reference_site import ReferenceSite
from panel.annotation import Annotation
from panel.drug_summary import DrugSummary
from panel.gene_panel import GenePanel
from panel.haplotype import Haplotype
from panel.panel import Panel
from panel.rs_id_info import RsIdInfo
from panel.variant import Variant

Json = Any


class JsonParser(object):
    # toplevel
    PANEL_NAME_KEY = "panelName"
    PANEL_VERSION_KEY = "panelVersion"
    GENES_KEY = "genes"
    IGNORED_VARIANTS_KEY = "ignoredVariants"

    # per gene
    GENE_NAME_KEY = "gene"
    CHROMOSOME_V37_KEY = "chromosomeV37"
    CHROMOSOME_V38_KEY = "chromosomeV38"
    WILD_TYPE_HAPLOTYPE_NAME_KEY = "wildTypeHaplotype"
    TRANSCRIPT_ID_KEY = "canonicalTranscript"
    RSID_INFOS_KEY = "variants"
    HAPLOTYPES_KEY = "haplotypes"
    DRUG_SUMMARIES_KEY = "drugs"
    REF_SEQ_DIFF_KEY = "refSeqDifferenceAnnotations"

    # per annotation
    ANNOTATION_V37_KEY = "annotationV37"
    ANNOTATION_V38_KEY = "annotationV38"

    # per rs id
    RS_ID_KEY = "rsid"
    REFERENCE_ALLELE_V37_KEY = "referenceAlleleV37"
    REFERENCE_ALLELE_V38_KEY = "referenceAlleleV38"
    POSITION_V37_KEY = "positionV37"
    POSITION_V38_KEY = "positionV38"

    # per drug
    DRUG_NAME_KEY = "name"
    PRESCRIPTION_URL_KEY = "urlPrescriptionInfo"

    # per haplotype
    HAPLOTYPE_NAME_KEY = "haplotypeName"
    HAPLOTYPE_FUNCTION_KEY = "function"
    HAPLOTYPE_VARIANTS_KEY = "haplotypeVariants"

    # per variant
    VARIANTS_RS_ID_KEY = "rsid"
    ALT_ALLELE_V37_KEY = "altAlleleV37"
    ALT_ALLELE_V38_KEY = "altAlleleV38"

    def get_panel(self, data: Json) -> Panel:
        name = str(data[self.PANEL_NAME_KEY])
        version = str(data[self.PANEL_VERSION_KEY])
        gene_panels = {self.get_gene_panel(gene_panel_json) for gene_panel_json in data[self.GENES_KEY]}
        if self.IGNORED_VARIANTS_KEY in data.keys():
            ignored_variants = {
                self.get_ignored_variant(variant_json) for variant_json in data[self.IGNORED_VARIANTS_KEY]
            }
        else:
            ignored_variants = set()
        return Panel(name, version, gene_panels, ignored_variants)

    def get_ignored_variant(self, data: Json) -> DualCall:
        chromosome_v37 = data[self.CHROMOSOME_V37_KEY]
        chromosome_v38 = data[self.CHROMOSOME_V38_KEY]
        position_v37 = data[self.POSITION_V37_KEY]
        position_v38 = data[self.POSITION_V38_KEY]
        reference_allele_v37 = data[self.REFERENCE_ALLELE_V37_KEY]
        reference_allele_v38 = data[self.REFERENCE_ALLELE_V38_KEY]
        alt_allele_v37 = data[self.ALT_ALLELE_V37_KEY]
        alt_allele_v38 = data[self.ALT_ALLELE_V38_KEY]

        reference_site_v37 = ReferenceSite(GeneCoordinate(chromosome_v37, int(position_v37)), reference_allele_v37)
        reference_site_v38 = ReferenceSite(GeneCoordinate(chromosome_v38, int(position_v38)), reference_allele_v38)
        call = DualCall(reference_site_v37, reference_site_v38, alt_allele_v37, alt_allele_v38)
        return call

    def get_gene_panel(self, data: Json) -> GenePanel:
        gene = str(data[self.GENE_NAME_KEY])
        chromosome_v37 = str(data[self.CHROMOSOME_V37_KEY])
        chromosome_v38 = str(data[self.CHROMOSOME_V38_KEY])
        wild_type_haplotype = str(data[self.WILD_TYPE_HAPLOTYPE_NAME_KEY])
        transcript_id: Optional[str]
        if self.TRANSCRIPT_ID_KEY in data.keys():
            transcript_id = str(data[self.TRANSCRIPT_ID_KEY])
        else:
            transcript_id = None
        rs_id_infos = self.get_rs_id_infos(
            data[self.RSID_INFOS_KEY],
            data[self.REF_SEQ_DIFF_KEY],
            chromosome_v37,
            chromosome_v38,
        )
        haplotypes = frozenset({self.get_haplotype(haplotype_json) for haplotype_json in data[self.HAPLOTYPES_KEY]})
        drugs = frozenset({self.get_drug_summary(drug_json) for drug_json in data[self.DRUG_SUMMARIES_KEY]})
        gene_panel = GenePanel(
            gene,
            wild_type_haplotype,
            transcript_id,
            haplotypes,
            rs_id_infos,
            drugs,
        )
        return gene_panel

    def get_rs_id_infos(
            self,
            rs_id_infos_data: Json,
            ref_seq_difference_annotation_data: Json,
            chromosome_v37: str,
            chromosome_v38: str,
    ) -> FrozenSet[RsIdInfo]:
        rs_id_to_ref_seq_difference_annotation = {
            str(annotation_json[self.RS_ID_KEY]): self.get_annotation(annotation_json)
            for annotation_json in ref_seq_difference_annotation_data
        }
        rs_id_infos = frozenset(
            {
                self.get_rs_id_info(
                    rs_id_info_json,
                    chromosome_v37,
                    chromosome_v38,
                    rs_id_to_ref_seq_difference_annotation,
                )
                for rs_id_info_json in rs_id_infos_data
            }
        )
        return rs_id_infos

    def get_rs_id_info(
            self,
            data: Json,
            chromosome_v37: str,
            chromosome_v38: str,
            rs_id_to_ref_seq_difference_annotation: Dict[str, Annotation],
    ) -> RsIdInfo:
        rs_id = str(data[self.RS_ID_KEY])
        reference_allele_v37 = str(data[self.REFERENCE_ALLELE_V37_KEY])
        reference_allele_v38 = str(data[self.REFERENCE_ALLELE_V38_KEY])
        start_coordinate_v37 = GeneCoordinate(chromosome_v37, int(data[self.POSITION_V37_KEY]))
        start_coordinate_v38 = GeneCoordinate(chromosome_v38, int(data[self.POSITION_V38_KEY]))
        annotation: Optional[Annotation]
        if reference_allele_v37 != reference_allele_v38 and rs_id in rs_id_to_ref_seq_difference_annotation.keys():
            # ignore annotation when alleles are identical, since we ignore additional fields in the panel json
            annotation = rs_id_to_ref_seq_difference_annotation[rs_id]
        else:
            annotation = None
        info = RsIdInfo(
            rs_id,
            ReferenceSite(start_coordinate_v37, reference_allele_v37),
            ReferenceSite(start_coordinate_v38, reference_allele_v38),
            annotation,
        )
        return info

    def get_haplotype(self, data: Json) -> Haplotype:
        name = str(data[self.HAPLOTYPE_NAME_KEY])
        function = str(data[self.HAPLOTYPE_FUNCTION_KEY])
        variants = frozenset({self.get_variant(variant_json) for variant_json in data[self.HAPLOTYPE_VARIANTS_KEY]})
        return Haplotype(name, function, variants)

    def get_drug_summary(self, data: Json) -> DrugSummary:
        name = str(data[self.DRUG_NAME_KEY])
        url_prescription_info = str(data[self.PRESCRIPTION_URL_KEY])
        return DrugSummary(name, url_prescription_info)

    def get_annotation(self, data: Json) -> Annotation:
        annotation_v37 = str(data[self.ANNOTATION_V37_KEY])
        annotation_v38 = str(data[self.ANNOTATION_V38_KEY])
        return Annotation(annotation_v37, annotation_v38)

    def get_variant(self, data: Json) -> Variant:
        rs_id = str(data[self.VARIANTS_RS_ID_KEY])
        variant_allele = str(data[self.ALT_ALLELE_V38_KEY])
        return Variant(rs_id, variant_allele)
