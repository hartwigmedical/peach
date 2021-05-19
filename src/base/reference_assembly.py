from enum import Enum, auto


class ReferenceAssembly(Enum):
    V37 = auto()
    V38 = auto()

    def opposite(self) -> "ReferenceAssembly":
        if self == ReferenceAssembly.V37:
            return ReferenceAssembly.V38
        elif self == ReferenceAssembly.V38:
            return ReferenceAssembly.V37
        else:
            raise NotImplementedError(f"Unrecognized reference assembly: {self}")

