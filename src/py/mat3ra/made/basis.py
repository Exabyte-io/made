import json
from typing import Dict, List, Optional

from mat3ra.code.constants import AtomicCoordinateUnits
from mat3ra.utils.mixins import RoundNumericValuesMixin
from pydantic import BaseModel

from .cell import Cell
from .utils import ArrayWithIds


class Basis(RoundNumericValuesMixin, BaseModel):
    elements: ArrayWithIds = ArrayWithIds(values=["Si"])
    coordinates: ArrayWithIds = ArrayWithIds(values=[0, 0, 0])
    units: str = AtomicCoordinateUnits.crystal
    cell: Optional[Cell] = None
    labels: Optional[ArrayWithIds] = ArrayWithIds(values=[])
    constraints: Optional[ArrayWithIds] = ArrayWithIds(values=[])

    @classmethod
    def from_dict(
        cls,
        elements: List[Dict],
        coordinates: List[Dict],
        units: str,
        labels: Optional[List[Dict]] = None,
        cell: Optional[Dict] = None,
        constraints: Optional[List[Dict]] = None,
    ) -> "Basis":
        return Basis(
            elements=ArrayWithIds.from_list_of_dicts(elements),
            coordinates=ArrayWithIds.from_list_of_dicts(coordinates),
            units=units,
            cell=Cell.from_nested_array(cell),
            labels=ArrayWithIds.from_list_of_dicts(labels) if labels else ArrayWithIds(values=[]),
            constraints=ArrayWithIds.from_list_of_dicts(constraints) if constraints else ArrayWithIds(values=[]),
        )

    def to_json(self, skip_rounding=False):
        json_value = {
            "elements": self.elements.to_json(),
            "coordinates": self.coordinates.to_json(skip_rounding=skip_rounding),
            "units": self.units,
            "cell": self.cell.to_json(skip_rounding=skip_rounding) if self.cell else None,
            "labels": self.labels.to_json(),
        }
        return json.loads(json.dumps(json_value))

    def clone(self):
        return Basis(
            elements=self.toJSON()["elements"],
            coordinates=self.toJSON()["coordinates"],
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
        self.coordinates = self.coordinates.map_array_in_place(self.cell.convert_point_to_cartesian)
        self.units = AtomicCoordinateUnits.cartesian

    def to_crystal(self):
        if self.is_in_crystal_units:
            return
        self.coordinates = self.coordinates.map_array_in_place(self.cell.convert_point_to_crystal)
        self.units = AtomicCoordinateUnits.crystal

    def add_atom(self, element="Si", coordinate=[0.5, 0.5, 0.5]):
        self.elements.add_item(element)
        self.coordinates.add_item(coordinate)

    def remove_atom_by_id(self, id=None):
        self.elements.remove_item(id)
        self.coordinates.remove_item(id)
        self.labels.remove_item(id)

    def filter_atoms_by_ids(self, ids):
        self.elements.filter_by_ids(ids)
        self.coordinates.filter_by_ids(ids)
        if self.labels is not None:
            self.labels.filter_by_ids(ids)
