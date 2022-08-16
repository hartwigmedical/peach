from typing import Optional, FrozenSet, Set, Dict

from util.gene_coordinate import GeneCoordinate
from util.reference_assembly import ReferenceAssembly
from util.reference_site import ReferenceSite
from panel.annotation import Annotation
from panel.rs_id_info import RsIdInfo


class GeneTranscriptSummary(object):
    def __init__(
            self,
            gene: str,
            transcript_id: Optional[str],
            rs_id_infos: FrozenSet[RsIdInfo],
    ) -> None:
        self.__assert_rs_id_infos_do_not_overlap(rs_id_infos, gene)
        self.__assert_rs_id_infos_match_chromosome(rs_id_infos, gene)

        rs_id_to_info: Dict[str, RsIdInfo] = {}
        for info in rs_id_infos:
            if info.rs_id in rs_id_to_info.keys():
                error_msg = f"Rs id '{info.rs_id}' present multiple times for gene '{gene}'"
                raise ValueError(error_msg)
            rs_id_to_info[info.rs_id] = info

        self.__gene = gene
        self.__transcript_id = transcript_id
        self.__rs_id_to_info = rs_id_to_info

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, GeneTranscriptSummary)
                and self.__gene == other.__gene
                and self.__transcript_id == other.__transcript_id
                and self.__rs_id_to_info == other.__rs_id_to_info
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.__gene,
                self.__transcript_id,
                self.get_rs_id_infos(),
            )
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"GeneTranscriptSummary("
            f"gene={self.__gene!r}, "
            f"transcript_id={self.__transcript_id!r}, "
            f"rs_id_infos={self.get_rs_id_infos()!r}, "
            f")"
        )

    @property
    def gene(self) -> str:
        return self.__gene

    @property
    def transcript_id(self) -> Optional[str]:
        return self.__transcript_id

    def get_rs_ids(self) -> Set[str]:
        return set(self.__rs_id_to_info.keys())

    def get_ref_seq_difference_annotation(self, rs_id: str) -> Optional[Annotation]:
        return self.__rs_id_to_info[rs_id].ref_seq_difference_annotation

    def get_reference_site(self, rs_id: str, reference_assembly: ReferenceAssembly) -> ReferenceSite:
        return self.__rs_id_to_info[rs_id].get_reference_site(reference_assembly)

    def get_rs_id_infos(self) -> FrozenSet[RsIdInfo]:
        return frozenset(self.__rs_id_to_info.values())

    @staticmethod
    def __assert_rs_id_infos_do_not_overlap(rs_id_infos: FrozenSet[RsIdInfo], gene: str) -> None:
        seen_v37_coordinates: Set[GeneCoordinate] = set()
        seen_v38_coordinates: Set[GeneCoordinate] = set()
        for info in rs_id_infos:
            if seen_v37_coordinates.intersection(info.reference_site_v37.get_covered_coordinates()):
                error_msg = (
                    f"Some rs id infos for gene '{gene}' have overlapping v37 coordinates. "
                    f"One of them has rs id '{info.rs_id}'."
                )
                raise ValueError(error_msg)
            seen_v37_coordinates = seen_v37_coordinates.union(info.reference_site_v37.get_covered_coordinates())

            if seen_v38_coordinates.intersection(info.reference_site_v38.get_covered_coordinates()):
                error_msg = (
                    f"Some rs id infos for gene '{gene}' have overlapping v38 coordinates. "
                    f"One of them has rs id '{info.rs_id}'."
                )
                raise ValueError(error_msg)
            seen_v38_coordinates = seen_v38_coordinates.union(info.reference_site_v38.get_covered_coordinates())

    @staticmethod
    def __assert_rs_id_infos_match_chromosome(rs_id_infos: FrozenSet[RsIdInfo], gene: str) -> None:
        v37_chromosomes = {info.reference_site_v37.start_coordinate.chromosome for info in rs_id_infos}
        v38_chromosomes = {info.reference_site_v38.start_coordinate.chromosome for info in rs_id_infos}
        if len(v37_chromosomes) > 1 or len(v38_chromosomes) > 1:
            error_msg = (
                f"Rs id infos for gene {gene} disagree on chromosome:\n"
                f"v37 chromosomes: {v37_chromosomes}\n"
                f"v38 chromosomes: {v38_chromosomes}"
            )
            raise ValueError(error_msg)
