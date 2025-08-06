from typing import Any, Type

from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseSingleBuilder

from .......modify import translate_to_z_level, wrap_to_unit_cell
from .......operations.core.unary import perturb
from ...... import MaterialWithBuildMetadata
from ..build_parameters import PerturbationBuildParameters
from ..configuration import PerturbationConfiguration


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
        self, material: MaterialWithBuildMetadata, configuration: PerturbationConfiguration
    ) -> MaterialWithBuildMetadata:
        perturbation_details = f"Perturbation: {configuration.perturbation_function_holder.function_str}"
        material.name = f"{material.name} ({perturbation_details})"
        return material
