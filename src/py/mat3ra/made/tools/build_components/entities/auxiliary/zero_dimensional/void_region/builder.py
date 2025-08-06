from typing import Type

from ......build_components import MaterialWithBuildMetadata
from ......modify import filter_by_condition_on_coordinates
from ..point_defect_site.builder import PointDefectSiteBuilder
from ..void_region.configuration import VoidRegionConfiguration


class VoidRegionBuilder(PointDefectSiteBuilder):
    """
    Builder class for creating a material with a void region based on a coordinate condition.
    Uses "Vac" as the element for all atoms that fall outside the specified coordinate condition
        and those positions will void elements in target crystal when merging.
    """

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
