class HaplotypeCall(object):
    def __init__(self, haplotype_name: str, count: int) -> None:
        if not 1 <= count <= 2:
            error_msg = f"Illegal haplotype count {count} for haplotype {haplotype_name}"
            raise SyntaxError(error_msg)

        self.__haplotype_name = haplotype_name
        self.__count = count

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, HaplotypeCall)
            and self.__haplotype_name == other.__haplotype_name
            and self.__count == other.__count
        )

    def __hash__(self) -> int:
        return hash((self.__haplotype_name, self.__count))

    def __repr__(self) -> str:  # pragma: no cover
        return f"HaplotypeCall({self.__haplotype_name!r}, {self.__count!r})"

    @property
    def haplotype_name(self) -> str:
        return self.__haplotype_name

    @property
    def count(self) -> int:
        return self.__count
