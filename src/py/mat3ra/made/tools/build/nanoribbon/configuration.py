from typing import Literal

from mat3ra.made.material import Material
from pydantic import BaseModel
from mat3ra.code.entity import InMemoryEntity


class NanoribbonConfiguration(BaseModel, InMemoryEntity):
    """
    Configuration for building a nanoribbon.

    Attributes:
        material (Material): The material to build the nanoribbon from.
        width (int): The width of the nanoribbon in number of unit cells.
        length (int): The length of the nanoribbon in number of unit cells.
        edge_type (Literal["armchair", "zigzag"]): The edge type of the nanoribbon.
    """

    material: Material
    width: int  # in number of unit cells
    length: int  # in number of unit cells
    vacuum_width: int = 3  # in number of unit cells
    vacuum_length: int = 0  # in number of unit cells
    edge_type: Literal["armchair", "zigzag"]

    class Config:
        arbitrary_types_allowed = True

    @property
    def _json(self):
        return {
            "material": self.material.to_json(),
            "width": self.width,
            "length": self.length,
            "vacuum_width": self.vacuum_width,
            "vacuum_length": self.vacuum_length,
            "edge_type": self.edge_type,
        }
