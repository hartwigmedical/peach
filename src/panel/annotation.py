from typing import NamedTuple

from util.reference_assembly import ReferenceAssembly


class Annotation(NamedTuple):
    annotation_v37: str
    annotation_v38: str

    def for_assembly(self, reference_assembly: ReferenceAssembly) -> str:
        if reference_assembly == ReferenceAssembly.V37:
            return self.annotation_v37
        elif reference_assembly == ReferenceAssembly.V38:
            return self.annotation_v38
        else:
            raise NotImplementedError(f"Unrecognized reference assembly: {reference_assembly}")
