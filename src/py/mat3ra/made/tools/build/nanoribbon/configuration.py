from mat3ra.made.material import Material
from mat3ra.made.tools.build.nanoribbon.enums import EdgeTypes

from ...build import BaseConfiguration


class NanoribbonConfiguration(BaseConfiguration):
    """
    Configuration for building a nanoribbon from a material.


    Attributes:
        material (Material): The material to build the nanoribbon from.
        width (int): The width of the nanoribbon in number of unit cells.
        length (int): The length of the nanoribbon in number of unit cells.
        vacuum_width (int): The width of the vacuum region in number of unit cells.
        vacuum_length (int): The length of the vacuum region in number of unit cells.
        edge_type (EdgeTypes): The type of edge to use for the nanoribbon, either zigzag or armchair.
    """

    material: Material
    width: int  # in number of unit cells
    length: int  # in number of unit cells
    vacuum_width: int = 3  # in number of unit cells
    vacuum_length: int = 0  # in number of unit cells
    edge_type: EdgeTypes = EdgeTypes.zigzag

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
