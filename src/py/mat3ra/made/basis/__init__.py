from typing import Any, Dict, List, Optional, Union

from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.material import BasisSchema
from mat3ra.esse.models.material import Units as UnitsEnum
from mat3ra.made.basis.coordinates import Coordinates
from mat3ra.made.cell import Cell
from mat3ra.made.utils import get_overlapping_coordinates
from pydantic import Field


class Basis(BasisSchema, InMemoryEntityPydantic):
    elements: ArrayWithIds
    coordinates: Coordinates
    cell: Cell = Field(Cell(), exclude=True)
    labels: ArrayWithIds = Field(ArrayWithIds.from_values([]))
    constraints: ArrayWithIds = Field(ArrayWithIds.from_values([]))

    def __convert_kwargs__(self, **kwargs: Any) -> Dict[str, Any]:
        if isinstance(kwargs.get("elements"), list):
            kwargs["elements"] = ArrayWithIds.from_list_of_dicts(kwargs["elements"])
        if isinstance(kwargs.get("coordinates"), list):
            kwargs["coordinates"] = Coordinates.from_list_of_dicts(kwargs["coordinates"])
        if isinstance(kwargs.get("labels"), list):
            kwargs["labels"] = ArrayWithIds.from_list_of_dicts(kwargs["labels"])
        if isinstance(kwargs.get("constraints"), list):
            kwargs["constraints"] = ArrayWithIds.from_list_of_dicts(kwargs["constraints"])
        if isinstance(kwargs.get("cell"), list):
            kwargs["cell"] = Cell.from_vectors_array(kwargs["cell"])
        return kwargs

    def __init__(self, *args: Any, **kwargs: Any):
        kwargs = self.__convert_kwargs__(**kwargs)
        super().__init__(*args, **kwargs)

    @classmethod
    def from_dict(
        cls,
        elements: List[Dict],
        coordinates: List[Dict],
        units: str,
        cell: List[List[float]],
        labels: Optional[List[Dict]] = ArrayWithIds.from_list_of_dicts([]),
        constraints: Optional[List[Dict]] = ArrayWithIds.from_list_of_dicts([]),
    ) -> "Basis":
        return Basis(
            elements=ArrayWithIds.from_list_of_dicts(elements),
            coordinates=Coordinates.from_list_of_dicts(coordinates),
            units=units,
            cell=Cell.from_vectors_array(cell),
            labels=ArrayWithIds.from_list_of_dicts(labels),
            constraints=ArrayWithIds.from_list_of_dicts(constraints),
        )

    @property
    def is_in_crystal_units(self):
        return self.units == UnitsEnum.crystal

    @property
    def is_in_cartesian_units(self):
        return self.units == UnitsEnum.cartesian

    def to_cartesian(self):
        if self.is_in_cartesian_units:
            return
        self.coordinates.map_array_in_place(self.cell.convert_point_to_cartesian)
        self.units = UnitsEnum.cartesian

    def to_crystal(self):
        if self.is_in_crystal_units:
            return
        self.coordinates.map_array_in_place(self.cell.convert_point_to_crystal)
        self.units = UnitsEnum.crystal

    def add_atom(
        self,
        element="Si",
        coordinate: Optional[List[float]] = None,
        use_cartesian_coordinates: bool = False,
        force: bool = False,
    ):
        """
        Add an atom to the basis.

        Before adding the atom at the specified coordinate, checks that no other atom is overlapping within a threshold.

        Args:
            element (str): Element symbol of the atom to be added.
            coordinate (List[float]): Coordinate of the atom to be added.
            use_cartesian_coordinates (bool): Whether the coordinate is in Cartesian units (or crystal by default).
            force (bool): Whether to force adding the atom even if it overlaps with another atom.
        """
        if coordinate is None:
            coordinate = [0, 0, 0]
        if use_cartesian_coordinates and self.is_in_crystal_units:
            coordinate = self.cell.convert_point_to_crystal(coordinate)
        if not use_cartesian_coordinates and self.is_in_cartesian_units:
            coordinate = self.cell.convert_point_to_cartesian(coordinate)
        cartesian_coordinates_for_overlap_check = [
            self.cell.convert_point_to_cartesian(coord) for coord in self.coordinates.values
        ]
        cartesian_coordinate_for_overlap_check = self.cell.convert_point_to_cartesian(coordinate)
        if get_overlapping_coordinates(
            cartesian_coordinate_for_overlap_check, cartesian_coordinates_for_overlap_check, threshold=0.01
        ):
            if force:
                print(f"Warning: Overlapping coordinates found for {coordinate}. Adding atom anyway.")
            else:
                print(f"Warning: Overlapping coordinates found for {coordinate}. Not adding atom.")
                return
        self.elements.add_item(element)
        self.coordinates.add_item(coordinate)

    def add_atoms_from_another_basis(self, other_basis: "Basis"):
        """
        Add atoms from another basis to this basis.

        Args:
            other_basis (Basis): The other basis to add atoms from.
        """

        self.elements.add_items(other_basis.elements.values)
        self.coordinates.add_items(other_basis.coordinates.values)
        self.labels.add_items(other_basis.labels.values)

    def remove_atom_by_id(self, id: int):
        self.elements.remove_item(id)
        self.coordinates.remove_item(id)
        self.labels.remove_item(id)

    def filter_atoms_by_ids(self, ids: Union[List[int], int], invert: bool = False) -> "Basis":
        self.elements.filter_by_ids(ids, invert)
        self.coordinates.filter_by_ids(ids, invert)
        self.labels.filter_by_ids(ids, invert)
        return self

    def filter_atoms_by_labels(self, labels: Union[List[str], str]) -> "Basis":
        self.labels.filter_by_values(labels)
        ids = self.labels.ids
        self.elements.filter_by_ids(ids)
        self.coordinates.filter_by_ids(ids)
        return self
