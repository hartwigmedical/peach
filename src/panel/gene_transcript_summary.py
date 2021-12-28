import itertools
from typing import Optional, FrozenSet, Set

from panel.rs_id_info import RsIdInfo


class GeneTranscriptSummary(object):
    def __init__(
            self,
            gene: str,
            transcript_id: Optional[str],
            rs_id_infos: FrozenSet[RsIdInfo],
    ) -> None:
        self.__assert_rs_ids_all_different(rs_id_infos)
        self.__assert_rs_id_infos_compatible(rs_id_infos)
        self.__assert_rs_id_infos_match_chromosome(rs_id_infos, gene)
    
        self.__gene = gene
        self.__transcript_id = transcript_id
        self.__rs_id_infos = rs_id_infos

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, GeneTranscriptSummary)
                and self.__gene == other.__gene
                and self.__transcript_id == other.__transcript_id
                and self.__rs_id_infos == other.__rs_id_infos
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.__gene,
                self.__transcript_id,
                self.__rs_id_infos,
            )
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"GeneTranscriptSummary("
            f"gene={self.__gene!r}, "
            f"transcript_id={self.__transcript_id!r}, "
            f"rs_id_infos={self.__rs_id_infos!r}, "
            f")"
        )

    @property
    def gene(self) -> str:
        return self.__gene

    @property
    def transcript_id(self) -> Optional[str]:
        return self.__transcript_id

    @property
    def rs_id_infos(self) -> FrozenSet[RsIdInfo]:
        return self.__rs_id_infos

    def get_rs_ids(self) -> Set[str]:
        return {rs_id_info.rs_id for rs_id_info in self.__rs_id_infos}

    @staticmethod
    def __assert_rs_id_infos_compatible(rs_id_infos: FrozenSet[RsIdInfo]) -> None:
        for left, right in itertools.combinations(rs_id_infos, 2):
            if not left.is_compatible(right):
                error_msg = f"Incompatible rs id infos in gene info. left: {left}, right: {right}"
                raise ValueError(error_msg)

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

    @staticmethod
    def __assert_rs_ids_all_different(rs_id_infos: FrozenSet[RsIdInfo]) -> None:
        rs_ids = [info.rs_id for info in rs_id_infos]
        if len(rs_ids) != len(set(rs_ids)):
            error_msg = (
                f"Not all rs ids are different: rs_ids={sorted(rs_ids)}"
            )
            raise ValueError(error_msg)
