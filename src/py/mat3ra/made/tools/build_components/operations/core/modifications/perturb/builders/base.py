from typing import Any, Type

from ... import MaterialWithBuildMetadata, BaseSingleBuilder, TypeConfiguration
from ..configuration import PerturbationConfiguration
from ..build_parameters import PerturbationBuildParameters
from ....modify import wrap_to_unit_cell, translate_to_z_level
from ....operations.core.unary import perturb


class PerturbationBuilder(BaseSingleBuilder):
    _ConfigurationType: Type[PerturbationConfiguration] = PerturbationConfiguration
    _BuildParametersType: Type[PerturbationBuildParameters] = PerturbationBuildParameters
    _DefaultBuildParameters: PerturbationBuildParameters = PerturbationBuildParameters()
    _PostProcessParametersType: Any = None

    def _generate(self, configuration: PerturbationConfiguration) -> MaterialWithBuildMetadata:
        new_material = translate_to_z_level(configuration.material, "center").clone()
        new_material = perturb(
            new_material,
            configuration.perturbation_function_holder,
            configuration.use_cartesian_coordinates,
        )
        new_material = wrap_to_unit_cell(new_material)
        return new_material

    def _update_material_name(
        self, material: MaterialWithBuildMetadata, configuration: TypeConfiguration
    ) -> MaterialWithBuildMetadata:
        perturbation_details = f"Perturbation: {configuration.perturbation_function_holder.function_str}"
        material.name = f"{material.name} ({perturbation_details})"
        return material
