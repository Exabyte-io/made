from typing import Dict, Type

from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.point.builders import PointDefectSiteBuilder
from mat3ra.made.tools.build.merge import MergeBuilder
from mat3ra.made.tools.build.void_region.configuration import VoidRegionConfiguration
from mat3ra.made.tools.modify import filter_by_condition_on_coordinates


class VoidRegionBuilder(PointDefectSiteBuilder):
    # we need to replace all atoms elements whose coordinates fall outside of the coordinate condition with "Vac" and remove the rest.
    _ConfigurationType: Type[VoidRegionConfiguration] = VoidRegionConfiguration

    def _generate(self, configuration: VoidRegionConfiguration) -> MaterialWithBuildMetadata:
        filtered_material = filter_by_condition_on_coordinates(
            configuration.crystal,
            configuration.coordinate_condition.condition,
            use_cartesian_coordinates=configuration.use_cartesian_coordinates,
            invert_selection=not configuration.invert_selection,
        )

        # Replace all elements with "Vac"
        filtered_material.basis.elements.values = ["Vac"] * len(filtered_material.basis.elements.values)

        return filtered_material
