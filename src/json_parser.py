from typing import Optional

from base.gene_coordinate import GeneCoordinate
from base.json_alias import Json
from base.reference_site import ReferenceSite
from config.annotation import Annotation
from config.drug_info import DrugInfo
from config.gene_info import GeneInfo
from config.haplotype import Haplotype
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant


class JsonParser(object):
    @classmethod
    def get_panel(cls, data: Json) -> Panel:
        name = str(data["panelName"])
        version = str(data["panelVersion"])
        gene_infos = frozenset({cls.get_gene_info(gene_info_json) for gene_info_json in data["genes"]})
        return Panel(name, version, gene_infos)

    @classmethod
    def get_gene_info(cls, data: Json) -> GeneInfo:
        gene = str(data["gene"])
        chromosome_v37 = str(data["chromosomeV37"])
        chromosome_v38 = str(data["chromosomeV38"])
        wild_type_haplotype = str(data["wildTypeHaplotype"])
        transcript_id: Optional[str]
        if "canonicalTranscript" in data.keys():
            transcript_id = str(data["canonicalTranscript"])
        else:
            transcript_id = None
        rs_id_infos = frozenset(
            {
                cls.get_rs_id_info(rs_id_info_json, chromosome_v37, chromosome_v38)
                for rs_id_info_json in data["variants"]
            }
        )
        haplotypes = frozenset({cls.get_haplotype(haplotype_json) for haplotype_json in data["haplotypes"]})
        drugs = frozenset({cls.get_drug_info(drug_json) for drug_json in data["drugs"]})
        rs_id_to_ref_seq_difference_annotation_v38 = {
            str(annotation_json["rsid"]): cls.get_annotation(annotation_json)
            for annotation_json in data["refSeqDifferenceAnnotations"]
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
        rs_id = str(data["rsid"])
        reference_allele_v37 = str(data["referenceAlleleV37"])
        reference_allele_v38 = str(data["referenceAlleleV38"])
        start_coordinate_v37 = GeneCoordinate(chromosome_v37, int(data["positionV37"]))
        start_coordinate_v38 = GeneCoordinate(chromosome_v38, int(data["positionV38"]))
        info = RsIdInfo(
            rs_id,
            ReferenceSite(start_coordinate_v37, reference_allele_v37),
            ReferenceSite(start_coordinate_v38, reference_allele_v38),
        )
        return info

    @classmethod
    def get_haplotype(cls, data: Json) -> Haplotype:
        name = str(data["haplotypeName"])
        function = str(data["function"])
        variants = frozenset({cls.get_variant(variant_json) for variant_json in data["haplotypeVariants"]})
        return Haplotype(name, function, variants)

    @classmethod
    def get_drug_info(cls, data: Json) -> DrugInfo:
        name = str(data["name"])
        url_prescription_info = str(data["urlPrescriptionInfo"])
        return DrugInfo(name, url_prescription_info)

    @classmethod
    def get_annotation(cls, data: Json) -> Annotation:
        annotation_v37 = str(data["annotationV37"])
        annotation_v38 = str(data["annotationV38"])
        return Annotation(annotation_v37, annotation_v38)

    @classmethod
    def get_variant(cls, data: Json) -> Variant:
        rs_id = str(data["rsid"])
        variant_allele = str(data["altAlleleV38"])
        return Variant(rs_id, variant_allele)
