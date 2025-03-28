from typing import Any, List, ClassVar, Dict

from mat3ra.code.constants import AtomicCoordinateUnits, Units
from mat3ra.code.entity import HasDescriptionHasMetadataNamedDefaultableInMemoryEntityPydantic
from mat3ra.esse.models.material import MaterialSchema, BasisSchema

from .basis import Basis
from .lattice import Lattice

defaultMaterialConfig = {
    "name": "Silicon FCC",
    "basis": {
        "elements": [
            {
                "id": 1,
                "value": "Si",
            },
            {
                "id": 2,
                "value": "Si",
            },
        ],
        "coordinates": [
            {
                "id": 1,
                "value": [0.0, 0.0, 0.0],
            },
            {
                "id": 2,
                "value": [0.25, 0.25, 0.25],
            },
        ],
        "units": AtomicCoordinateUnits.crystal,
    },
    "lattice": {
        "type": "FCC",
        "a": 3.867,
        "b": 3.867,
        "c": 3.867,
        "alpha": 60,
        "beta": 60,
        "gamma": 60,
        "units": {
            "length": Units.angstrom,
            "angle": Units.degree,
        },
    },
}


# TODO: replace `-Pydantic` with actual class in the next PR
class Material(MaterialSchema, HasDescriptionHasMetadataNamedDefaultableInMemoryEntityPydantic):
    __default_config__: ClassVar[Dict[str, Any]] = defaultMaterialConfig

    def model_post_init(self, __context: Any) -> None:
        if not self.name and self.formula:
            self.name = self.formula

    @property
    def coordinates_array(self) -> List[List[float]]:
        return self.BasisCls.coordinates.values

    @property
    def BasisCls(self) -> Basis:
        config = self.model_dump()["basis"]
        config["cell"] = config.get("cell", self.LatticeCls.vector_arrays)
        return Basis.from_dict(**config)

    @property
    def LatticeCls(self) -> Lattice:
        return Lattice(**self.model_dump()["lattice"])

    def to_cartesian(self) -> None:
        new_basis = Basis(**self.basis.model_dump())  # convert from BasisSchema → Basis
        new_basis.to_cartesian()
        self.basis = BasisSchema(**new_basis.model_dump())

    def to_crystal(self) -> None:
        new_basis = self.basis.model_copy()
        new_basis.to_crystal()
        self.basis = BasisSchema(**new_basis.model_dump())

    def set_coordinates(self, coordinates: List[List[float]]) -> None:
        new_basis = self.basis.model_copy()
        new_basis.coordinates.values = coordinates
        self.basis = BasisSchema(**self.basis.model_dump())

    def set_new_lattice_vectors(
        self, lattice_vector1: List[float], lattice_vector2: List[float], lattice_vector3: List[float]
    ) -> None:
        lattice = Lattice.from_vectors_array([lattice_vector1, lattice_vector2, lattice_vector3])
        original_is_in_crystal = self.basis.is_in_crystal_units
        self.to_cartesian()
        self.lattice = lattice.to_schema()
        if original_is_in_crystal:
            self.to_crystal()

    def add_atom(self, element: str, coordinate: List[float], use_cartesian_coordinates: bool = False) -> None:
        new_basis = self.basis.model_copy()
        new_basis.add_atom(element, coordinate, use_cartesian_coordinates)
        self.basis = BasisSchema(**new_basis.model_dump())
