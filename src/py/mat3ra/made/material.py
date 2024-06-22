from typing import Any, Dict, List, Union

from mat3ra.code.constants import AtomicCoordinateUnits, Units
from mat3ra.code.entity import HasDescriptionHasMetadataNamedDefaultableInMemoryEntity
from mat3ra.esse.models.material import MaterialSchema

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

MaterialSchemaJSON = Dict[str, Union[MaterialSchema, Any]]


class Material(HasDescriptionHasMetadataNamedDefaultableInMemoryEntity):
    jsonSchema: MaterialSchemaJSON
    default_config = defaultMaterialConfig

    def __init__(self, config: Any) -> None:
        super().__init__(config)
        self.name = super().name or self.formula

    def to_json(self, exclude: List[str] = []) -> MaterialSchemaJSON:
        return {**super().to_json()}

    @property
    def coordinates_array(self) -> List[List[float]]:
        return self.basis.coordinates.values

    @property
    def basis(self) -> Basis:
        return Basis.from_dict(**self.get_prop("basis"))

    @basis.setter
    def basis(self, basis: Basis) -> None:
        self.set_prop("basis", basis.to_json())

    @property
    def lattice(self) -> Lattice:
        return Lattice(**self.get_prop("lattice"))

    @lattice.setter
    def lattice(self, lattice: Lattice) -> None:
        self.set_prop("lattice", lattice.to_json())
