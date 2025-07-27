from typing import Tuple, Optional

from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.miller_indices import (
    MillerIndicesSchema,
)
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.termination import (
    TerminationSchema,
)


class Termination(TerminationSchema, InMemoryEntityPydantic):
    def __str__(self):
        return f"{self.chemical_elements}_{self.space_group_symmetry_label}"

    def __repr__(self):
        return self.__str__()

    def __init__(self, chemical_elements: str, space_group_symmetry_label: str):
        super().__init__(chemical_elements=chemical_elements, space_group_symmetry_label=space_group_symmetry_label)

    def __eq__(self, other):
        if not isinstance(other, Termination):
            return NotImplemented
        return (
            self.chemical_elements == other.chemical_elements
            and self.space_group_symmetry_label == other.space_group_symmetry_label
        )

    @property
    def formula(self):
        return self.chemical_elements

    @classmethod
    def from_string(cls, termination: str):
        chemical_elements = termination.split("_")[0]
        space_group_symmetry_label = "_".join(termination.split("_")[1:])
        return cls(chemical_elements=chemical_elements, space_group_symmetry_label=space_group_symmetry_label)


class TerminationHolder(InMemoryEntityPydantic):
    termination_with_vacuum: Termination
    termination_without_vacuum: Optional[Termination]
    shift_with_vacuum: float
    shift_without_vacuum: Optional[float]


class MillerIndices(MillerIndicesSchema, InMemoryEntityPydantic):
    def to_tuple(self) -> Tuple[int, int, int]:
        return (self.root[0], self.root[1], self.root[2])

    @staticmethod
    def format_index(index: int) -> str:
        return f"{index}" + chr(0x0302) if index < 0 else f"{index}"

    def __str__(self) -> str:
        return f"({''.join([self.format_index(i) for i in self.root])})"
