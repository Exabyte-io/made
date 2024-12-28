import json
from typing import Dict, List, Optional, Union

from mat3ra.code.constants import AtomicCoordinateUnits
from mat3ra.utils.mixins import RoundNumericValuesMixin
from pydantic import BaseModel

from .cell import Cell
from .utils import ArrayWithIds, get_overlapping_coordinates


class Basis(RoundNumericValuesMixin, BaseModel):
    elements: ArrayWithIds = ArrayWithIds(values=["Si"])
    coordinates: ArrayWithIds = ArrayWithIds(values=[0, 0, 0])
    units: str = AtomicCoordinateUnits.crystal
    cell: Cell = Cell()
    labels: Optional[ArrayWithIds] = ArrayWithIds(values=[])
    constraints: Optional[ArrayWithIds] = ArrayWithIds(values=[])

    @classmethod
    def from_dict(
        cls,
        elements: List[Dict],
        coordinates: List[Dict],
        units: str,
        labels: Optional[List[Dict]] = None,
        cell: Optional[List[List[float]]] = None,
        constraints: Optional[List[Dict]] = None,
    ) -> "Basis":
        return Basis(
            elements=ArrayWithIds.from_list_of_dicts(elements),
            coordinates=ArrayWithIds.from_list_of_dicts(coordinates),
            units=units,
            cell=Cell.from_vectors_array(cell),
            labels=ArrayWithIds.from_list_of_dicts(labels) if labels else ArrayWithIds(values=[]),
            constraints=ArrayWithIds.from_list_of_dicts(constraints) if constraints else ArrayWithIds(values=[]),
        )

    def to_json(self, skip_rounding=False):
        json_value = {
            "elements": self.elements.to_json(),
            "coordinates": self.coordinates.to_json(skip_rounding=skip_rounding),
            "units": self.units,
            "labels": self.labels.to_json(),
        }
        return json.loads(json.dumps(json_value))

    def clone(self):
        return Basis(
            elements=self.elements,
            coordinates=self.coordinates,
            units=self.units,
            cell=self.cell,
            isEmpty=False,
            labels=self.labels,
        )

    @property
    def is_in_crystal_units(self):
        return self.units == AtomicCoordinateUnits.crystal

    @property
    def is_in_cartesian_units(self):
        return self.units == AtomicCoordinateUnits.cartesian

    def to_cartesian(self):
        if self.is_in_cartesian_units:
            return
        self.coordinates.map_array_in_place(self.cell.convert_point_to_cartesian)
        self.units = AtomicCoordinateUnits.cartesian

    def to_crystal(self):
        if self.is_in_crystal_units:
            return
        self.coordinates.map_array_in_place(self.cell.convert_point_to_crystal)
        self.units = AtomicCoordinateUnits.crystal

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

    def remove_atom_by_id(self, id: int):
        self.elements.remove_item(id)
        self.coordinates.remove_item(id)
        if self.labels is not None:
            self.labels.remove_item(id)

    def filter_atoms_by_ids(self, ids: Union[List[int], int]) -> "Basis":
        self.elements.filter_by_ids(ids)
        self.coordinates.filter_by_ids(ids)
        if self.labels is not None:
            self.labels.filter_by_ids(ids)
        return self

    def filter_atoms_by_labels(self, labels: Union[List[str], str]) -> "Basis":
        if self.labels is None:
            return self
        self.labels.filter_by_values(labels)
        ids = self.labels.ids
        self.elements.filter_by_ids(ids)
        self.coordinates.filter_by_ids(ids)
        return self
