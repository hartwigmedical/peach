from typing import Optional

from base.reference_site import ReferenceSite
from base.reference_assembly import ReferenceAssembly
from panel.annotation import Annotation


class RsIdInfo(object):
    def __init__(
            self,
            rs_id: str,
            reference_site_v37: ReferenceSite,
            reference_site_v38: ReferenceSite,
            ref_seq_difference_annotation: Optional[Annotation],
    ) -> None:
        if reference_site_v37.allele == reference_site_v38.allele and ref_seq_difference_annotation is not None:
            error_msg = (
                f"For rs id '{rs_id}' a 'reference sequence difference annotation' is provided "
                f"even though the reference sequences are not different in this location: "
                f"allele_v37={reference_site_v37.allele}, "
                f"allele_v38={reference_site_v38.allele}, "
                f"annotation={ref_seq_difference_annotation}"
            )
            raise ValueError(error_msg)
        if reference_site_v37.allele != reference_site_v38.allele and ref_seq_difference_annotation is None:
            error_msg = (
                f"For rs id '{rs_id}' a 'reference sequence difference annotation' has not been provided "
                f"even though the reference sequences are different in this location: "
                f"allele_v37={reference_site_v37.allele}, "
                f"allele_v38={reference_site_v38.allele}, "
                f"annotation={ref_seq_difference_annotation}"
            )
            raise ValueError(error_msg)
        self.__rs_id = rs_id
        self.__reference_site_v37 = reference_site_v37
        self.__reference_site_v38 = reference_site_v38
        self.__ref_seq_difference_annotation = ref_seq_difference_annotation

    @property
    def rs_id(self) -> str:
        return self.__rs_id

    @property
    def reference_site_v37(self) -> ReferenceSite:
        return self.__reference_site_v37

    @property
    def reference_site_v38(self) -> ReferenceSite:
        return self.__reference_site_v38

    @property
    def ref_seq_difference_annotation(self) -> Optional[Annotation]:
        return self.__ref_seq_difference_annotation

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, RsIdInfo)
                and self.__rs_id == other.__rs_id
                and self.__reference_site_v37 == other.__reference_site_v37
                and self.__reference_site_v38 == other.__reference_site_v38
                and self.__ref_seq_difference_annotation == other.__ref_seq_difference_annotation
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.__rs_id,
                self.__reference_site_v37,
                self.__reference_site_v38,
                self.__ref_seq_difference_annotation,
            )
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"RsIdInfo("
            f"rs_id={self.__rs_id!r}, "
            f"reference_site_v37={self.__reference_site_v37!r}, "
            f"reference_site_v38={self.__reference_site_v38!r}, "
            f"annotation={self.__ref_seq_difference_annotation!r}, "
            f")"
        )

    def get_reference_site(self, reference_assembly: ReferenceAssembly) -> ReferenceSite:
        if reference_assembly == ReferenceAssembly.V37:
            return self.reference_site_v37
        elif reference_assembly == ReferenceAssembly.V38:
            return self.reference_site_v38
        else:
            error_msg = "Unrecognized reference assembly version"
            raise NotImplementedError(error_msg)
