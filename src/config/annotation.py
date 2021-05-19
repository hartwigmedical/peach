from typing import NamedTuple

from base.json_alias import Json


class Annotation(NamedTuple):
    annotation_v37: str
    annotation_v38: str

    @classmethod
    def from_json(cls, data: Json) -> "Annotation":
        annotation_v37 = str(data["annotationV37"])
        annotation_v38 = str(data["annotationV38"])
        return Annotation(annotation_v37, annotation_v38)
