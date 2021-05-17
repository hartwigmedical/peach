from typing import NamedTuple

from base.json_alias import Json


class Annotation(NamedTuple):
    annotation_v38: str

    @classmethod
    def from_json(cls, data: Json) -> "Annotation":
        return Annotation(str(data["annotationV38"]))
