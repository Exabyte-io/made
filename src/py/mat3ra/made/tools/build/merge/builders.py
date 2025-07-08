from typing import Any, Optional, List

from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseSingleBuilder, BaseBuilderParameters
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration
from mat3ra.made.tools.build.vacuum.builders import VacuumBuilder
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.made.tools.operations.core.binary import merge


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

    def _merge_add(self, materials: List[Material], parameters: MergeBuilderParameters) -> Material:
        return merge(
            materials=materials,
            merge_method=MergeMethodsEnum.ADD,
            material_name=parameters.material_name or materials[0].name,
            distance_tolerance=parameters.distance_tolerance,
            merge_dangerously=parameters.merge_dangerously,
        )

    def _merge_replace(self, materials: List[Material], parameters: MergeBuilderParameters) -> Material:
        return merge(
            materials=materials,
            merge_method=MergeMethodsEnum.REPLACE,
            material_name=parameters.material_name or materials[0].name,
            distance_tolerance=parameters.distance_tolerance,
            merge_dangerously=parameters.merge_dangerously,
        )

    def _merge_yield(self, materials: List[Material], parameters: MergeBuilderParameters) -> Material:
        return merge(
            materials=materials,
            merge_method=MergeMethodsEnum.YIELD,
            material_name=parameters.material_name or materials[0].name,
            distance_tolerance=parameters.distance_tolerance,
            merge_dangerously=parameters.merge_dangerously,
        )

    def _generate(self, configuration: MergeConfiguration) -> Material:
        materials = []
        for component in configuration.merge_components:
            material = self._configuration_to_material(component)
            materials.append(material)

        if len(materials) == 1:
            return materials[0]

        parameters = self.build_parameters or self._DefaultBuildParameters
        merge_method = getattr(configuration, "merge_method", None)
        if merge_method is not None:
            if merge_method == MergeMethodsEnum.ADD:
                return self._merge_add(materials, parameters)
            elif merge_method == MergeMethodsEnum.REPLACE:
                return self._merge_replace(materials, parameters)
            elif merge_method == MergeMethodsEnum.YIELD:
                return self._merge_yield(materials, parameters)
            else:
                raise ValueError(f"Unknown merge method: {merge_method}")
        else:
            return self._merge_add(materials, parameters)
