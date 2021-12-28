from enum import Enum, auto


class VcfCallFilter(Enum):
    # For calls with respect to the vcf reference assembly
    PASS = auto()
    NO_CALL = auto()


class DualCallFilter(Enum):
    # For dual calls
    PASS = auto()
    NO_CALL = auto()
    UNKNOWN = auto()
    INFERRED_PASS = auto()
