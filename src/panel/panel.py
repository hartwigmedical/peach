from typing import Set, FrozenSet, Dict, Optional, Union

from util.gene_coordinate import GeneCoordinate
from util.reference_assembly import ReferenceAssembly
from util.reference_site import ReferenceSite
from calls.single_call import SingleCall
from calls.vcf_call import VcfCall
from panel.gene_panel import GenePanel
from panel.variant import Variant


class Panel(object):
    def __init__(self, name: str, version: str, gene_panels: FrozenSet[GenePanel]) -> None:
        gene_to_gene_panel: Dict[str, GenePanel] = {}
        for gene_panel in gene_panels:
            if gene_panel.gene in gene_to_gene_panel.keys():
                error_msg = f"The panel '{name}' contains the gene '{gene_panel.gene}' more than once."
                raise ValueError(error_msg)
            gene_to_gene_panel[gene_panel.gene] = gene_panel

        rs_id_to_gene: Dict[str, str] = {}
        for gene_panel in gene_panels:
            for rs_id in gene_panel.get_rs_ids():
                if rs_id in rs_id_to_gene.keys():
                    error_msg = f"The rs id {rs_id} is included for multiple genes."
                    raise ValueError(error_msg)
                rs_id_to_gene[rs_id] = gene_panel.gene

        self.__name = name
        self.__version = version
        self.__gene_to_gene_panel = gene_to_gene_panel

        self.__rs_id_to_gene = rs_id_to_gene
        self.__v37_reference_site_to_rs_id = self.__get_reference_site_to_rs_id(
            gene_panels, ReferenceAssembly.V37, name
        )
        self.__v38_reference_site_to_rs_id = self.__get_reference_site_to_rs_id(
            gene_panels, ReferenceAssembly.V38, name
        )

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, Panel)
                and self.__name == other.__name
                and self.__version == other.__version
                and self.__gene_to_gene_panel == other.__gene_to_gene_panel
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"Panel("
            f"name={self.__name!r}, "
            f"version={self.__version!r}, "
            f"gene_panels={frozenset(self.__gene_to_gene_panel.values())!r}, "
            f")"
        )

    @property
    def name(self) -> str:
        return self.__name

    @property
    def version(self) -> str:
        return self.__version

    def get_drugs(self, gene: str) -> Set[str]:
        return self.__gene_to_gene_panel[gene].get_drug_names()

    def get_drug_prescription_url(self, gene: str, drug_name: str) -> str:
        return self.__gene_to_gene_panel[gene].get_prescription_url(drug_name)

    def get_reference_site(self, rs_id: str, reference_assembly: ReferenceAssembly) -> ReferenceSite:
        return self.__gene_to_gene_panel[self.get_gene_for_rs_id(rs_id)].get_reference_site(rs_id, reference_assembly)

    def get_genes(self) -> Set[str]:
        return set(self.__gene_to_gene_panel.keys())

    def get_rs_ids_for_gene(self, gene: str) -> Set[str]:
        return self.__gene_to_gene_panel[gene].get_rs_ids()

    def get_wild_type_haplotype_name(self, gene: str) -> str:
        return self.__gene_to_gene_panel[gene].wild_type_haplotype_name

    def get_non_wild_type_haplotype_names(self, gene: str) -> Set[str]:
        return self.__gene_to_gene_panel[gene].get_non_wild_type_haplotype_names()

    def get_variants_for_haplotype(self, gene: str, haplotype_name: str) -> Set[Variant]:
        return self.__gene_to_gene_panel[gene].get_variants(haplotype_name)

    def get_haplotype_function(self, gene: str, haplotype_name: str) -> str:
        return self.__gene_to_gene_panel[gene].get_haplotype_function(haplotype_name)

    def get_transcript_ids(self) -> Set[str]:
        transcript_ids = {
            gene_panel.transcript_id for gene_panel in self.__gene_to_gene_panel.values()
            if gene_panel.transcript_id is not None
        }
        return transcript_ids

    def is_relevant_to_panel(self, call: VcfCall, reference_assembly: ReferenceAssembly) -> bool:
        reference_site_to_rs_id = self.__get_matching_reference_site_to_rs_id(reference_assembly)
        covered_coordinates_call = call.reference_site.get_covered_coordinates()

        for reference_site, rs_id in reference_site_to_rs_id.items():
            if reference_site.get_covered_coordinates().intersection(covered_coordinates_call) or rs_id in call.rs_ids:
                return True
        return False

    def get_relevant_panel_rs_ids(self, call: VcfCall, reference_assembly: ReferenceAssembly) -> Set[str]:
        reference_site_to_rs_id = self.__get_matching_reference_site_to_rs_id(reference_assembly)
        covered_coordinates_call = call.reference_site.get_covered_coordinates()

        relevant_rs_ids = set()
        for reference_site, rs_id in reference_site_to_rs_id.items():
            if reference_site.get_covered_coordinates().intersection(covered_coordinates_call) or rs_id in call.rs_ids:
                relevant_rs_ids.add(rs_id)

        return relevant_rs_ids

    def get_perfectly_matching_rs_id(
            self, call: Union[VcfCall, SingleCall], reference_assembly: ReferenceAssembly
    ) -> Optional[str]:
        matching_rs_id: Optional[str]
        if self.__contains_rs_id_with_reference_site(call.reference_site, reference_assembly):
            matching_rs_id = self.__get_rs_id_with_reference_site(call.reference_site, reference_assembly)
        else:
            if any(self.__contains_rs_id(rs_id) for rs_id in call.rs_ids):
                error_msg = f"Rs id is in panel, but the call reference site does not match the panel reference site."
                raise ValueError(error_msg)
            matching_rs_id = None
        return matching_rs_id

    def has_ref_seq_difference_annotation(self, rs_id: str) -> bool:
        annotation = self.__gene_to_gene_panel[self.get_gene_for_rs_id(rs_id)].get_ref_seq_difference_annotation(rs_id)
        return annotation is not None

    def get_ref_seq_difference_annotation(self, rs_id: str, reference_assembly: ReferenceAssembly) -> str:
        annotation = self.__gene_to_gene_panel[self.get_gene_for_rs_id(rs_id)].get_ref_seq_difference_annotation(rs_id)
        if annotation is None:
            raise ValueError(f"No ref seq difference annotation for rs_id: rs_id={rs_id}")
        return annotation.for_assembly(reference_assembly)

    def get_gene_for_rs_id(self, rs_id: str) -> str:
        return self.__rs_id_to_gene[rs_id]

    def is_empty(self) -> bool:
        return not self.__gene_to_gene_panel.keys()

    def get_id(self) -> str:
        return f"{self.__name}_v{self.__version}"

    def __contains_rs_id_with_reference_site(
            self, reference_site: ReferenceSite, reference_assembly: ReferenceAssembly
    ) -> bool:
        return reference_site in self.__get_matching_reference_site_to_rs_id(reference_assembly)

    def __get_rs_id_with_reference_site(
            self, reference_site: ReferenceSite, reference_assembly: ReferenceAssembly
    ) -> str:
        return self.__get_matching_reference_site_to_rs_id(reference_assembly)[reference_site]

    def __contains_rs_id(self, rs_id: str) -> bool:
        return rs_id in self.__rs_id_to_gene.keys()

    def __get_matching_reference_site_to_rs_id(self, reference_assembly: ReferenceAssembly) -> Dict[ReferenceSite, str]:
        if reference_assembly == ReferenceAssembly.V37:
            return self.__v37_reference_site_to_rs_id
        elif reference_assembly == ReferenceAssembly.V38:
            return self.__v38_reference_site_to_rs_id
        else:
            raise ValueError(f"Unrecognized reference assembly: '{reference_assembly}'")

    @staticmethod
    def __get_reference_site_to_rs_id(
            gene_panels: FrozenSet[GenePanel], reference_assembly: ReferenceAssembly, name: str
    ) -> Dict[ReferenceSite, str]:
        seen_covered_coordinates: Set[GeneCoordinate] = set()
        reference_site_to_rs_id: Dict[ReferenceSite, str] = {}
        for gene_panel in gene_panels:
            for rs_id in gene_panel.get_rs_ids():
                reference_site = gene_panel.get_reference_site(rs_id, reference_assembly)
                if seen_covered_coordinates.intersection(reference_site.get_covered_coordinates()):
                    error_msg = (
                        f"Some rs ids for panel '{name}' have overlapping {reference_assembly.name} coordinates. "
                        f"One of them is rs id '{rs_id}' for gene '{gene_panel.gene}'."
                    )
                    raise ValueError(error_msg)
                seen_covered_coordinates = seen_covered_coordinates.union(reference_site.get_covered_coordinates())
                reference_site_to_rs_id[reference_site] = rs_id
        return reference_site_to_rs_id
