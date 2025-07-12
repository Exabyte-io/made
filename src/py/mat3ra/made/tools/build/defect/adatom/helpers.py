from typing import List, Optional

from mat3ra.esse.models.materials_category_components.entities.auxiliary.zero_dimensional.point_defect_site import (
    AtomSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.adatom import AdatomMaterialAnalyzer
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from .builders import AdatomDefectBuilder
from .configuration import (
    AdatomDefectConfiguration,
)
from ..point.builders import AtomAtCoordinateBuilder
from ..point.configuration import PointDefectSite
from ..slab.helpers import recreate_slab_with_fractional_layers
from ...defect.enums import AdatomPlacementMethodEnum


def create_adatom_defect(
    slab: Material,
    position_on_surface: List[float],
    distance_z: float = 1.0,
    placement_method: AdatomPlacementMethodEnum = AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE,
    element: Optional[str] = None,
) -> Material:
    """
    Create an adatom defect using the new AdatomDefectConfiguration and AdatomDefectBuilder.

    Args:
        slab: The slab material.
        position_on_surface: Position on the surface [x, y].
        element: Chemical element for the adatom.
        distance_z: Distance above the surface in Angstroms.
        placement_method: Method for placing the adatom.

    Returns:
        Material: The slab with adatom defect.
    """
    coordinate = [*position_on_surface, 0]

    slab_analyzer = SlabMaterialAnalyzer(material=slab, coordinate=coordinate)
    added_slab = recreate_slab_with_fractional_layers(slab, 1)

    adatom_analyzer = AdatomMaterialAnalyzer(
        material=slab,
        coordinate=coordinate,
        placement_method=placement_method,
        distance_z=distance_z,
    )

    resolved_coordinate = adatom_analyzer.coordinate_in_added_component
    slab_without_vacuum_configuration = adatom_analyzer.slab_without_vacuum_configuration

    adatom_point_defect_config = PointDefectSite(
        crystal=added_slab,
        element=AtomSchema(chemical_element=element),
        coordinate=resolved_coordinate,
    )
    isolated_defect = AtomAtCoordinateBuilder().get_material(adatom_point_defect_config)

    vacuum = slab_analyzer.get_slab_vacuum_configuration()

    configuration = AdatomDefectConfiguration(
        stack_components=[slab_without_vacuum_configuration, isolated_defect, vacuum],
    )

    builder = AdatomDefectBuilder()
    return builder.get_material(configuration)
