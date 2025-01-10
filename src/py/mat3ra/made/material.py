from typing import Any, Dict, List, Union, Optional

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
        config = self.get_prop("basis")
        config["cell"] = config.get("cell", self.lattice.vector_arrays)
        return Basis.from_dict(**config)

    @basis.setter
    def basis(self, basis: Basis) -> None:
        self.set_prop("basis", basis.to_json())

    @property
    def lattice(self) -> Lattice:
        return Lattice(**self.get_prop("lattice"))

    @lattice.setter
    def lattice(self, lattice: Lattice) -> None:
        self.set_prop("lattice", lattice.to_json())

    def to_cartesian(self) -> None:
        new_basis = self.basis.copy()
        new_basis.to_cartesian()
        self.basis = new_basis

    def to_crystal(self) -> None:
        new_basis = self.basis.copy()
        new_basis.to_crystal()
        self.basis = new_basis

    def set_coordinates(self, coordinates: List[List[float]]) -> None:
        new_basis = self.basis.copy()
        new_basis.coordinates.values = coordinates
        self.basis = new_basis

    def set_new_lattice_vectors(
        self, lattice_vector1: List[float], lattice_vector2: List[float], lattice_vector3: List[float]
    ) -> None:
        lattice = Lattice.from_vectors_array([lattice_vector1, lattice_vector2, lattice_vector3])
        original_is_in_crystal = self.basis.is_in_crystal_units
        self.to_cartesian()
        self.lattice = lattice
        if original_is_in_crystal:
            self.to_crystal()

    def add_atom(self, element: str, coordinate: List[float], use_cartesian_coordinates=False) -> None:
        new_basis = self.basis.copy()
        new_basis.add_atom(element, coordinate, use_cartesian_coordinates)
        self.basis = new_basis


    @classmethod
    def create_empty(
        cls,
        a: float = 1.0,
        b: Optional[float] = None,
        c: Optional[float] = None,
        alpha: float = 90.0,
        beta: float = 90.0,
        gamma: float = 90.0,
        lattice_type: str = "CUB",
        name: str = "New Material"
    ) -> "Material":
        """
        Create an empty material with specified lattice parameters.

        Args:
            a (float): First lattice constant. Defaults to 1.0.
            b (float, optional): Second lattice constant. Defaults to value of 'a'.
            c (float, optional): Third lattice constant. Defaults to value of 'a'.
            alpha (float): First lattice angle in degrees. Defaults to 90.0.
            beta (float): Second lattice angle in degrees. Defaults to 90.0.
            gamma (float): Third lattice angle in degrees. Defaults to 90.0.
            lattice_type (str): Type of lattice. Defaults to "CUB".
            name (str): Name of the material. Defaults to "New Material".

        Returns:
            Material: A new empty Material instance with specified lattice
        """
        b = b if b is not None else a
        c = c if c is not None else a

        basis_config = {
            "elements": [],
            "coordinates": [],
            "units": AtomicCoordinateUnits.cartesian,
        }

        lattice_config = {
            "type": lattice_type,
            "a": a,
            "b": b,
            "c": c,
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
            "units": {
                "length": Units.angstrom,
                "angle": Units.degree,
            }
        }

        config = {
            "name": name,
            "basis": basis_config,
            "lattice": lattice_config
        }

        return cls(config)
