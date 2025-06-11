from typing import Tuple, Optional

from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.material.reusable.slab.miller_indices import MillerIndicesSchema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.termination import (
    TerminationSchema,
)
from pydantic import BaseModel


class Termination(TerminationSchema, InMemoryEntityPydantic):
    def __str__(self):
        return f"{self.chemical_elements}_{self.space_group_symmetry_label}"

    def __repr__(self):
        return self.__str__()

    def __init__(self, chemical_elements: str, space_group_symmetry_label: str):
        super().__init__(chemical_elements=chemical_elements, space_group_symmetry_label=space_group_symmetry_label)

    def __eq__(self, other):
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


class TerminationHolder(BaseModel):
    termination_with_vacuum: Termination
    termination_without_vacuum: Optional[Termination]
    shift_with_vacuum: float
    shift_without_vacuum: Optional[float]


class MillerIndices(MillerIndicesSchema, InMemoryEntityPydantic):
    def to_tuple(self) -> Tuple[int, int, int]:
        return (self.root[0], self.root[1], self.root[2])
