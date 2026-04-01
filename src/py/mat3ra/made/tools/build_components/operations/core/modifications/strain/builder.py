from typing import Type

from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseBuilderParameters, BaseSingleBuilder

from ..... import MaterialWithBuildMetadata
from ......operations.core.unary import strain
from .configuration import StrainConfiguration


class StrainBuilder(BaseSingleBuilder):
    _ConfigurationType: Type[StrainConfiguration] = StrainConfiguration
    _BuildParametersType: Type[BaseBuilderParameters] = BaseBuilderParameters
    _DefaultBuildParameters: BaseBuilderParameters = BaseBuilderParameters()

    def _generate(self, configuration: StrainConfiguration) -> MaterialWithBuildMetadata:
        strained_material = strain(configuration.material, configuration.strain_matrix)
        return MaterialWithBuildMetadata.create(strained_material.to_dict())
