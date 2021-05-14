from enum import Enum, auto


class SimpleCallFilter(Enum):
    # For v37 calls
    PASS = auto()
    NO_CALL = auto()


class FullCallFilter(Enum):
    # For full calls
    PASS = auto()
    NO_CALL = auto()
    UNKNOWN = auto()
    INFERRED_PASS = auto()
