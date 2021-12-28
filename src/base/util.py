from typing import TypeVar

S = TypeVar("S")
T = TypeVar("T")


def replace_file_extension_of_path(path: str, new_file_extension: str) -> str:
    assert "." in path, f"Path does not contain point to signal file extension: path={path}"
    split_path = path.split(".")
    new_path = ".".join(split_path[0:-1]) + "." + new_file_extension
    return new_path


def strip_prefix(string: str, prefix: str) -> str:
    if string.startswith(prefix):
        return string[len(prefix):]
    else:
        error_msg = (
            f"String does not start with expected prefix: string={string}, prefix={prefix}"
        )
        raise ValueError(error_msg)
