from typing import Optional, Any

from base.gene_coordinate import GeneCoordinate
from base.reference_site import ReferenceSite
from config.annotation import Annotation
from config.drug_info import DrugInfo
from config.gene_info import GeneInfo
from config.haplotype import Haplotype
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant

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

    @classmethod
    def get_panel(cls, data: Json) -> Panel:
        name = str(data[cls.PANEL_NAME])
        version = str(data[cls.PANEL_VERSION])
        gene_infos = frozenset({cls.get_gene_info(gene_info_json) for gene_info_json in data[cls.GENES_KEY]})
        return Panel(name, version, gene_infos)

    @classmethod
    def get_gene_info(cls, data: Json) -> GeneInfo:
        gene = str(data[cls.GENE_NAME])
        chromosome_v37 = str(data[cls.CHROMOSOME_V37])
        chromosome_v38 = str(data[cls.CHROMOSOME_V38])
        wild_type_haplotype = str(data[cls.WILD_TYPE_HAPLOTYPE_NAME])
        transcript_id: Optional[str]
        if cls.TRANSCRIPT_ID in data.keys():
            transcript_id = str(data[cls.TRANSCRIPT_ID])
        else:
            transcript_id = None
        rs_id_infos = frozenset(
            {
                cls.get_rs_id_info(rs_id_info_json, chromosome_v37, chromosome_v38)
                for rs_id_info_json in data[cls.RSID_INFOS_KEY]
            }
        )
        haplotypes = frozenset({cls.get_haplotype(haplotype_json) for haplotype_json in data[cls.HAPLOTYPES_KEY]})
        drugs = frozenset({cls.get_drug_info(drug_json) for drug_json in data[cls.DRUG_INFOS_KEY]})
        rs_id_to_ref_seq_difference_annotation_v38 = {
            str(annotation_json[cls.RS_ID_KEY]): cls.get_annotation(annotation_json)
            for annotation_json in data[cls.REF_SEQ_DIFF_KEY]
        }
        gene_info = GeneInfo(
            gene,
            wild_type_haplotype,
            transcript_id,
            haplotypes,
            rs_id_infos,
            drugs,
            rs_id_to_ref_seq_difference_annotation_v38,
        )
        return gene_info

    @classmethod
    def get_rs_id_info(cls, data: Json, chromosome_v37: str, chromosome_v38: str) -> RsIdInfo:
        rs_id = str(data[cls.RS_ID_KEY])
        reference_allele_v37 = str(data[cls.REFERENCE_ALLELE_V37])
        reference_allele_v38 = str(data[cls.REFERENCE_ALLELE_V38])
        start_coordinate_v37 = GeneCoordinate(chromosome_v37, int(data[cls.POSITION_V37]))
        start_coordinate_v38 = GeneCoordinate(chromosome_v38, int(data[cls.POSITION_V38]))
        info = RsIdInfo(
            rs_id,
            ReferenceSite(start_coordinate_v37, reference_allele_v37),
            ReferenceSite(start_coordinate_v38, reference_allele_v38),
        )
        return info

    @classmethod
    def get_haplotype(cls, data: Json) -> Haplotype:
        name = str(data[cls.HAPLOTYPE_NAME])
        function = str(data[cls.HAPLOTYPE_FUNCTION])
        variants = frozenset({cls.get_variant(variant_json) for variant_json in data[cls.HAPLOTYPE_VARIANTS_KEY]})
        return Haplotype(name, function, variants)

    @classmethod
    def get_drug_info(cls, data: Json) -> DrugInfo:
        name = str(data[cls.DRUG_NAME])
        url_prescription_info = str(data[cls.PRESCRIPTION_URL])
        return DrugInfo(name, url_prescription_info)

    @classmethod
    def get_annotation(cls, data: Json) -> Annotation:
        annotation_v37 = str(data[cls.ANNOTATION_V37])
        annotation_v38 = str(data[cls.ANNOTATION_V38])
        return Annotation(annotation_v37, annotation_v38)

    @classmethod
    def get_variant(cls, data: Json) -> Variant:
        rs_id = str(data[cls.VARIANTS_RS_ID])
        variant_allele = str(data[cls.ALT_ALLELE_V38])
        return Variant(rs_id, variant_allele)
