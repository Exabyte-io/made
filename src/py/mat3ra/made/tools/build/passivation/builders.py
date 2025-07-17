from typing import List, Dict, Type
from mat3ra.made.material import Material
from mat3ra.made.tools.bonds import BondDirections
from pydantic import BaseModel
import numpy as np

from .. import MaterialWithBuildMetadata
from ..defect.point.builders import AtomAtCoordinateConfiguration, AtomAtCoordinateBuilder
from ..merge import MergeBuilder
from ...analyze.material import MaterialWithCrystalSites
from ...enums import SurfaceTypes
from ...analyze.other import get_surface_atom_indices
from ...modify import translate_to_z_level
from ...build import BaseBuilder
from .configuration import PassivationConfiguration


class PassivationBuilder(MergeBuilder):

    @property
    def merge_component_types_conversion_map(self) -> Dict[Type, Type]:
        return {
            AtomAtCoordinateConfiguration: AtomAtCoordinateBuilder,
        }

    def _update_material_name(
        self, material: BaseBuilder._GeneratedItemType, configuration: PassivationConfiguration
    ) -> MaterialWithBuildMetadata:
        material_name = configuration.material.name
        material.name = f"{material_name}, {configuration.passivant}-passivated"
        return material

    # def create_passivated_material(self, configuration: BaseBuilder._ConfigurationType) -> Material:
    #     material = configuration.material.clone()
    #     # material.to_cartesian()
    #     material = translate_to_z_level(material, "center")
    #     # material.to_crystal()
    #     return material
    #
    # def _add_passivant_atoms(
    #     self, material: Material, coordinates: list, passivant: str, use_cartesian_coordinates=False
    # ) -> Material:
    #     """
    #     Add passivant atoms to the provided coordinates in the material.
    #
    #     Args:
    #         material (Material): The material object to add passivant atoms to.
    #         coordinates (list): The coordinates to add passivant atoms to.
    #         passivant (str): The chemical symbol of the passivating atom (e.g., 'H').
    #         use_cartesian_coordinates (bool): Whether the coordinates are in Cartesian units (or crystal by default).
    #
    #     Returns:
    #         Material: The material object with passivation atoms added.
    #     """
    #     for coord in coordinates:
    #         material.add_atom(passivant, coord, use_cartesian_coordinates)
    #     return material
