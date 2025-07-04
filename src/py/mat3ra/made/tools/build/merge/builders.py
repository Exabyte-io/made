from typing import Any, Optional
from typing import TypeVar

import numpy as np
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseSingleBuilder, BaseBuilderParameters
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration
from mat3ra.made.tools.build.vacuum.builders import VacuumBuilder
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.made.tools.operations.core.binary import merge_materials
from mat3ra.made.tools.operations.reusable.unary import strain_to_match_lattice


class MergeBuilderParameters(BaseBuilderParameters):
    """
    Parameters for merging materials.

    Args:
        material_name: Optional name for the merged material.
        distance_tolerance: Tolerance for resolving overlapping coordinates.
        merge_dangerously: If True, allows merging even if lattices are different.
    """

    material_name: Optional[str] = "New Material"
    distance_tolerance: float = 0.1
    merge_dangerously: bool = False


class MergeBuilder(BaseSingleBuilder):
    """
    Builder class for merging materials based on parameters.
    Behaves like StackBuilder but performs merge operations.
    """

    _BuildParametersType = MergeBuilderParameters()
    _DefaultBuildParameters = MergeBuilderParameters()

    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, Material):
            return configuration_or_material

        if isinstance(configuration_or_material, VacuumConfiguration):
            builder = VacuumBuilder()
            return builder.get_material(configuration_or_material)
        raise ValueError(f"Unknown configuration type: {type(configuration_or_material)}")

    def _generate(self, configuration: MergeConfiguration) -> Material:
        materials = []
        for component in configuration.merge_components:
            material = self._configuration_to_material(component)
            materials.append(material)

        if len(materials) == 1:
            return materials[0]

        parameters = self.build_parameters or self._DefaultBuildParameters

        if parameters.merge_dangerously and len(materials) > 1:
            strained_materials = [materials[0]]
            for material_to_strain in materials[1:]:
                strained_material = strain_to_match_lattice(material_to_strain, materials[0])
                strained_materials.append(strained_material)

            materials = strained_materials

        merged_material = merge_materials(
            materials=materials,
            material_name=parameters.material_name,
            distance_tolerance=parameters.distance_tolerance,
            merge_dangerously=parameters.merge_dangerously,
        )
        return merged_material

