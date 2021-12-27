import itertools
from typing import Set, FrozenSet, Dict

from base.reference_assembly import ReferenceAssembly
from base.reference_site import ReferenceSite
from calls.vcf_call import VcfCall
from panel.gene_info import GeneInfo
from panel.haplotype import Haplotype
from panel.rs_id_info import RsIdInfo


class Panel(object):
    def __init__(self, name: str, version: str, gene_infos: FrozenSet[GeneInfo]) -> None:
        self.__assert_rs_ids_all_different(gene_infos)
        self.__assert_all_rs_id_infos_compatible(gene_infos)

        gene_to_gene_info: Dict[str, GeneInfo] = {}
        for gene_info in gene_infos:
            if gene_info.gene in gene_to_gene_info.keys():
                error_msg = f"The panel '{name}' contains the gene '{gene_info.gene}' more than once."
                raise ValueError(error_msg)
            gene_to_gene_info[gene_info.gene] = gene_info

        self.__name = name
        self.__version = version
        self.__gene_to_gene_info = gene_to_gene_info

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Panel)
            and self.__name == other.__name
            and self.__version == other.__version
            and self.__gene_to_gene_info == other.__gene_to_gene_info
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"Panel("
            f"name={self.__name!r}, "
            f"version={self.__version!r}, "
            f"gene_infos={frozenset(self.__gene_to_gene_info.values())!r}, "
            f")"
        )

    @property
    def name(self) -> str:
        return self.__name

    @property
    def version(self) -> str:
        return self.__version

    def get_drugs(self, gene: str) -> Set[str]:
        return {drug_info.name for drug_info in self.__gene_to_gene_info[gene].drugs}

    def get_drug_prescription_url(self, gene: str, drug_name: str) -> str:
        return self.__gene_to_gene_info[gene].get_prescription_url(drug_name)

    def get_reference_site(self, rs_id: str, reference_assembly: ReferenceAssembly) -> ReferenceSite:
        return self.__get_rs_id_info(rs_id).get_reference_site(reference_assembly)

    def get_genes(self) -> Set[str]:
        return set(self.__gene_to_gene_info.keys())

    def get_rs_ids_for_gene(self, gene: str) -> Set[str]:
        return self.__gene_to_gene_info[gene].get_rs_ids()

    def get_wild_type_haplotype_name(self, gene: str) -> str:
        return self.__gene_to_gene_info[gene].wild_type_haplotype_name

    def get_haplotypes(self, gene: str) -> Set[Haplotype]:
        return set(self.__gene_to_gene_info[gene].haplotypes)

    def get_transcript_ids(self) -> Set[str]:
        transcript_ids = {
            gene_info.transcript_id for gene_info in self.__gene_to_gene_info.values()
            if gene_info.transcript_id is not None
        }
        return transcript_ids

    def is_relevant_to_panel(self, call: VcfCall, reference_assembly: ReferenceAssembly) -> bool:
        covered_coordinates_call = call.reference_site.get_covered_coordinates()
        for info in self.__get_rs_id_infos():
            rs_id_covered_coordinates = info.get_reference_site(reference_assembly).get_covered_coordinates()
            if rs_id_covered_coordinates.intersection(covered_coordinates_call) or info.rs_id in call.rs_ids:
                return True
        return False

    def get_relevant_rs_ids(self, call: VcfCall, reference_assembly: ReferenceAssembly) -> Set[str]:
        relevant_rs_ids = set()
        covered_coordinates_call = call.reference_site.get_covered_coordinates()
        for info in self.__get_rs_id_infos():
            rs_id_covered_coordinates = info.get_reference_site(reference_assembly).get_covered_coordinates()
            if rs_id_covered_coordinates.intersection(covered_coordinates_call) or info.rs_id in call.rs_ids:
                relevant_rs_ids.add(info.rs_id)
        return relevant_rs_ids

    def contains_rs_id_with_reference_site(
            self, reference_site: ReferenceSite, reference_assembly: ReferenceAssembly
    ) -> bool:
        for info in self.__get_rs_id_infos():
            if info.get_reference_site(reference_assembly) == reference_site:
                return True
        return False

    def get_rs_id_with_reference_site(
            self, reference_site: ReferenceSite, reference_assembly: ReferenceAssembly
    ) -> str:
        matching_rs_ids = set()
        for info in self.__get_rs_id_infos():
            if info.get_reference_site(reference_assembly) == reference_site:
                matching_rs_ids.add(info.rs_id)

        if matching_rs_ids and len(matching_rs_ids) == 1:
            return matching_rs_ids.pop()
        elif not matching_rs_ids:
            raise ValueError("No rs id infos match position and ref allele")
        else:
            raise ValueError("Multiple rs id infos match position and ref allele. This should be impossible.")

    def has_ref_seq_difference_annotation(self, rs_id: str) -> bool:
        return self.__get_rs_id_info(rs_id).ref_seq_difference_annotation is not None

    def get_ref_seq_difference_annotation(self, rs_id: str, reference_assembly: ReferenceAssembly) -> str:
        annotation = self.__get_rs_id_info(rs_id).ref_seq_difference_annotation
        if annotation is None:
            raise ValueError(f"No ref seq difference annotation for rs_id: rs_id={rs_id}")
        return annotation.for_assembly(reference_assembly)

    def get_gene_for_rs_id(self, rs_id: str) -> str:
        matching_genes = []
        for gene_info in self.__gene_to_gene_info.values():
            for rs_id_info in gene_info.rs_id_infos:
                if rs_id_info.rs_id == rs_id:
                    matching_genes.append(gene_info.gene)
        if len(matching_genes) == 1:
            return matching_genes[0]
        else:
            raise ValueError(f"Not exactly one gene for rs_id: rs_id={rs_id}, genes={sorted(matching_genes)}")

    def contains_rs_id(self, rs_id: str) -> bool:
        return rs_id in self.__get_rs_ids()

    def is_empty(self) -> bool:
        return not self.__gene_to_gene_info.keys()

    def get_haplotype_function(self, gene: str, haplotype_name: str) -> str:
        return self.__gene_to_gene_info[gene].get_haplotype_function(haplotype_name)

    def get_id(self) -> str:
        return f"{self.__name}_v{self.__version}"

    def __get_rs_id_info(self, rs_id: str) -> RsIdInfo:
        matching_rs_id_infos = [rs_id_info for rs_id_info in self.__get_rs_id_infos() if rs_id_info.rs_id == rs_id]
        if len(matching_rs_id_infos) == 1:
            return matching_rs_id_infos[0]
        else:
            raise ValueError(f"Not exactly one matching rs id info in panel: rs_id={rs_id}")

    def __get_rs_ids(self) -> Set[str]:
        return {info.rs_id for info in self.__get_rs_id_infos()}

    def __get_rs_id_infos(self) -> Set[RsIdInfo]:
        return {rs_id_info for gene_info in self.__gene_to_gene_info.values() for rs_id_info in gene_info.rs_id_infos}

    @staticmethod
    def __assert_all_rs_id_infos_compatible(gene_infos: FrozenSet[GeneInfo]) -> None:
        for left_gene_info, right_gene_info in itertools.combinations(gene_infos, 2):
            for left_info in left_gene_info.rs_id_infos:
                for right_info in right_gene_info.rs_id_infos:
                    if not left_info.is_compatible(right_info):
                        error_msg = f"Incompatible rs id infos in panel. left: {left_info}, right: {right_info}"
                        raise ValueError(error_msg)

    @staticmethod
    def __assert_rs_ids_all_different(gene_infos: FrozenSet[GeneInfo]) -> None:
        rs_ids = [rs_id_info.rs_id for gene_info in gene_infos for rs_id_info in gene_info.rs_id_infos]
        if len(rs_ids) != len(set(rs_ids)):
            error_msg = (
                f"Not all rs ids are different: rs_ids={sorted(rs_ids)}"
            )
            raise ValueError(error_msg)
