from typing import Tuple

from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.miller_indices import (
    MillerIndicesSchema,
)


class MillerIndices(MillerIndicesSchema, InMemoryEntityPydantic):
    def to_tuple(self) -> Tuple[int, int, int]:
        return (self.root[0], self.root[1], self.root[2])

    @staticmethod
    def format_index(index: int) -> str:
        return f"{index}" + chr(0x0302) if index < 0 else f"{index}"

    def __str__(self) -> str:
        return f"({''.join([self.format_index(i) for i in self.root])})"
