from typing import Union

from mat3ra.made.material import Material
from .configuration import (
    PointDefectConfigurationLegacy,
    AdatomSlabPointDefectConfiguration,
)
from ...build import BaseBuilder


class DefectBuilder(BaseBuilder):
    def create_isolated_defect(self, defect_configuration: PointDefectConfigurationLegacy) -> Material:
        raise NotImplementedError


class DefectPairBuilder(DefectBuilder):
    def create_defect_pair(
        self,
        primary_defect_configuration: Union[PointDefectConfigurationLegacy, AdatomSlabPointDefectConfiguration],
        secondary_defect_configuration: Union[PointDefectConfigurationLegacy, AdatomSlabPointDefectConfiguration],
    ) -> Material:
        """
        Create a pair of point defects in the material.

        Args:
            primary_defect_configuration: The configuration of the first defect.
            secondary_defect_configuration: The configuration of the second defect.

        Returns:
            Material: The material with both defects added.
        """
        primary_material = self.create_isolated_defect(primary_defect_configuration)
        # Remove metadata to allow for independent defect creation
        if hasattr(primary_defect_configuration.crystal.metadata, "build"):
            primary_material.metadata["build"] = primary_defect_configuration.crystal.metadata["build"]
        primary_material.name = primary_defect_configuration.crystal.name
        secondary_defect_configuration.crystal = primary_material
        secondary_material = self.create_isolated_defect(secondary_defect_configuration)

        return secondary_material
