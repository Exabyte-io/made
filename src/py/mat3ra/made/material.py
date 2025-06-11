from typing import Any, List

from mat3ra.code.constants import AtomicCoordinateUnits, Units
from mat3ra.code.entity import HasDescriptionHasMetadataNamedDefaultableInMemoryEntityPydantic
from mat3ra.esse.models.material import MaterialSchema

from .basis import Basis
from .lattice import Lattice

defaultMaterialConfig = {
    "name": "Silicon FCC",
    "basis": {
        "elements": [
            {
                "id": 0,
                "value": "Si",
            },
            {
                "id": 1,
                "value": "Si",
            },
        ],
        "coordinates": [
            {
                "id": 0,
                "value": [0.0, 0.0, 0.0],
            },
            {
                "id": 1,
                "value": [0.25, 0.25, 0.25],
            },
        ],
        "units": AtomicCoordinateUnits.crystal,
        "labels": [],
        "constraints": [],
    },
    "lattice": {
        "a": 3.867,
        "b": 3.867,
        "c": 3.867,
        "alpha": 60.0,
        "beta": 60.0,
        "gamma": 60.0,
        "units": {
            "length": Units.angstrom,
            "angle": Units.degree,
        },
        "type": "FCC",
    },
}


# TODO: replace `-Pydantic` with actual class in the next PR
class Material(MaterialSchema, HasDescriptionHasMetadataNamedDefaultableInMemoryEntityPydantic):
    __default_config__ = defaultMaterialConfig
    __schema__ = MaterialSchema

    basis: Basis
    lattice: Lattice

    def model_post_init(self, __context: Any) -> None:
        if not self.name and self.formula:
            self.name: str = self.formula
        self.basis.cell = self.lattice.vectors

    @property
    def coordinates_array(self) -> List[List[float]]:
        return self.basis.coordinates.values

    def to_cartesian(self) -> None:
        self.basis.to_cartesian()

    def to_crystal(self) -> None:
        self.basis.to_crystal()

    def set_coordinates(self, coordinates: List[List[float]]) -> None:
        self.basis.coordinates.values = coordinates

    def set_new_lattice_vectors(
        self, lattice_vector1: List[float], lattice_vector2: List[float], lattice_vector3: List[float]
    ) -> None:
        original_is_in_crystal_units = self.basis.is_in_crystal_units
        self.to_cartesian()
        self.lattice = Lattice.from_vectors_array([lattice_vector1, lattice_vector2, lattice_vector3])
        self.basis.cell = self.lattice.vectors
        if original_is_in_crystal_units:
            self.to_crystal()

    def set_lattice(self, lattice: Lattice) -> None:
        self.set_new_lattice_vectors(*lattice.vector_arrays)

    def add_atom(self, element: str, coordinate: List[float], use_cartesian_coordinates: bool = False) -> None:
        self.basis.add_atom(element, coordinate, use_cartesian_coordinates)
