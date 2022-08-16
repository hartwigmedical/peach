from typing import TypeVar

S = TypeVar("S")
T = TypeVar("T")


def strip_prefix(string: str, prefix: str) -> str:
    if string.startswith(prefix):
        return string[len(prefix):]
    else:
        error_msg = (
            f"String does not start with expected prefix: string={string}, prefix={prefix}"
        )
        raise ValueError(error_msg)
