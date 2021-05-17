import inspect
from pathlib import Path


def __placeholder() -> None:  # pragma: no cover
    pass


def get_test_resources_dir() -> Path:
    source_file = inspect.getsourcefile(__placeholder)
    if source_file is None:  # pragma: no cover
        raise ValueError("Could not get test resource path.")
    return Path(source_file).parent


def get_test_resource(file_name: str) -> Path:
    return get_test_resources_dir() / file_name
