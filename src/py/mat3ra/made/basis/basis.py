import json
from typing import Dict

from mat3ra.code.constants import AtomicCoordinateUnits
from mat3ra.utils.mixins import RoundNumericValuesMixin
from pydantic import BaseModel

from ..cell.cell import Cell
from ..utils import ArrayWithIds


class Basis(RoundNumericValuesMixin, BaseModel):
    elements: ArrayWithIds = ArrayWithIds(array=["Si"])
    coordinates: ArrayWithIds = ArrayWithIds(array=[0, 0, 0])
    units: str = AtomicCoordinateUnits.crystal
    cell: Cell = Cell([1, 0, 0], [0, 1, 0], [0, 0, 1])
    # TODO: isolate labels to a separate class
    labels: ArrayWithIds = ArrayWithIds(array=[])

    @staticmethod
    def from_dict(config: Dict) -> "Basis":
        return Basis(
            elements=config.get("elements", []),
            coordinates=config.get("coordinates", []),
            units=config.get("units", AtomicCoordinateUnits.crystal),
            cell=Cell.from_dict(config.get("cell", {})),
            labels=config.get("labels", []),
        )

    def to_json(self, skip_rounding=False):
        json_value = {
            "elements": self.elements.to_json(),
            "coordinates": self.coordinates.to_json(skip_rounding=skip_rounding),
            "units": self.units,
            "cell": self.cell.to_json(skip_rounding=skip_rounding),
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
