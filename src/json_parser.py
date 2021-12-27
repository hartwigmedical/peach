from typing import Optional, Any, Dict, FrozenSet

from base.gene_coordinate import GeneCoordinate
from base.reference_site import ReferenceSite
from panel.annotation import Annotation
from panel.drug_info import DrugInfo
from panel.gene_info import GeneInfo
from panel.haplotype import Haplotype
from panel.panel import Panel
from panel.rs_id_info import RsIdInfo
from panel.variant import Variant

Json = Any


class JsonParser(object):
    # toplevel
    PANEL_NAME = "panelName"
    PANEL_VERSION = "panelVersion"
    GENES_KEY = "genes"

    # per gene
    GENE_NAME = "gene"
    CHROMOSOME_V37 = "chromosomeV37"
    CHROMOSOME_V38 = "chromosomeV38"
    WILD_TYPE_HAPLOTYPE_NAME = "wildTypeHaplotype"
    TRANSCRIPT_ID = "canonicalTranscript"
    RSID_INFOS_KEY = "variants"
    HAPLOTYPES_KEY = "haplotypes"
    DRUG_INFOS_KEY = "drugs"
    REF_SEQ_DIFF_KEY = "refSeqDifferenceAnnotations"

    # per annotation
    ANNOTATION_V37 = "annotationV37"
    ANNOTATION_V38 = "annotationV38"

    # per rs id
    RS_ID_KEY = "rsid"
    REFERENCE_ALLELE_V37 = "referenceAlleleV37"
    REFERENCE_ALLELE_V38 = "referenceAlleleV38"
    POSITION_V37 = "positionV37"
    POSITION_V38 = "positionV38"

    # per drug
    DRUG_NAME = "name"
    PRESCRIPTION_URL = "urlPrescriptionInfo"

    # per haplotype
    HAPLOTYPE_NAME = "haplotypeName"
    HAPLOTYPE_FUNCTION = "function"
    HAPLOTYPE_VARIANTS_KEY = "haplotypeVariants"

    # per variant
    VARIANTS_RS_ID = "rsid"
    ALT_ALLELE_V38 = "altAlleleV38"

    def get_panel(self, data: Json) -> Panel:
        name = str(data[self.PANEL_NAME])
        version = str(data[self.PANEL_VERSION])
        gene_infos = frozenset({self.get_gene_info(gene_info_json) for gene_info_json in data[self.GENES_KEY]})
        return Panel(name, version, gene_infos)

    def get_gene_info(self, data: Json) -> GeneInfo:
        gene = str(data[self.GENE_NAME])
        chromosome_v37 = str(data[self.CHROMOSOME_V37])
        chromosome_v38 = str(data[self.CHROMOSOME_V38])
        wild_type_haplotype = str(data[self.WILD_TYPE_HAPLOTYPE_NAME])
        transcript_id: Optional[str]
        if self.TRANSCRIPT_ID in data.keys():
            transcript_id = str(data[self.TRANSCRIPT_ID])
        else:
            transcript_id = None
        rs_id_infos = self.get_rs_id_infos(
            data[self.RSID_INFOS_KEY],
            data[self.REF_SEQ_DIFF_KEY],
            chromosome_v37,
            chromosome_v38,
        )
        haplotypes = frozenset({self.get_haplotype(haplotype_json) for haplotype_json in data[self.HAPLOTYPES_KEY]})
        drugs = frozenset({self.get_drug_info(drug_json) for drug_json in data[self.DRUG_INFOS_KEY]})
        gene_info = GeneInfo(
            gene,
            wild_type_haplotype,
            transcript_id,
            haplotypes,
            rs_id_infos,
            drugs,
        )
        return gene_info

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
        reference_allele_v37 = str(data[self.REFERENCE_ALLELE_V37])
        reference_allele_v38 = str(data[self.REFERENCE_ALLELE_V38])
        start_coordinate_v37 = GeneCoordinate(chromosome_v37, int(data[self.POSITION_V37]))
        start_coordinate_v38 = GeneCoordinate(chromosome_v38, int(data[self.POSITION_V38]))
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
        name = str(data[self.HAPLOTYPE_NAME])
        function = str(data[self.HAPLOTYPE_FUNCTION])
        variants = frozenset({self.get_variant(variant_json) for variant_json in data[self.HAPLOTYPE_VARIANTS_KEY]})
        return Haplotype(name, function, variants)

    def get_drug_info(self, data: Json) -> DrugInfo:
        name = str(data[self.DRUG_NAME])
        url_prescription_info = str(data[self.PRESCRIPTION_URL])
        return DrugInfo(name, url_prescription_info)

    def get_annotation(self, data: Json) -> Annotation:
        annotation_v37 = str(data[self.ANNOTATION_V37])
        annotation_v38 = str(data[self.ANNOTATION_V38])
        return Annotation(annotation_v37, annotation_v38)

    def get_variant(self, data: Json) -> Variant:
        rs_id = str(data[self.VARIANTS_RS_ID])
        variant_allele = str(data[self.ALT_ALLELE_V38])
        return Variant(rs_id, variant_allele)
