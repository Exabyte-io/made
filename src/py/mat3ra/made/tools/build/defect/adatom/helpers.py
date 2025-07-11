from typing import List, Optional

from mat3ra.esse.models.materials_category_components.entities.auxiliary.zero_dimensional.point_defect_site import (
    AtomSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from .builders import AdatomDefectBuilder
from .configuration import (
    AdatomDefectConfiguration,
)
from ..point.builders import PointDefectSiteBuilder
from ..point.configuration import PointDefectSite
from ..slab.helpers import recreate_slab_with_fractional_layers
from ...defect.enums import AdatomPlacementMethodEnum
from ....analyze.crystal_site import CrystalSiteAnalyzer
from ....analyze.other import get_atomic_coordinates_extremum


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
    distance_z_crystal = slab.basis.cell.convert_point_to_crystal([0, 0, distance_z])[2]
    coordinate = [*position_on_surface, distance_z_crystal]

    slab_analyzer = SlabMaterialAnalyzer(material=slab)

    slab_without_vacuum = slab_analyzer.get_slab_configuration_with_no_vacuum()
    added_slab = recreate_slab_with_fractional_layers(slab, 1)

    if placement_method == AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE:
        crystal_site_analyzer = CrystalSiteAnalyzer(material=added_slab, coordinate=coordinate)
        resolved_coordinate = crystal_site_analyzer.closest_site_coordinate
    elif placement_method == AdatomPlacementMethodEnum.EQUIDISTANT:
        max_z = get_atomic_coordinates_extremum(slab, "max", "z")
        coordinate[2] = max_z + distance_z_crystal
        crystal_site_analyzer = CrystalSiteAnalyzer(material=slab, coordinate=coordinate)
        resolved_coordinate = crystal_site_analyzer.equidistant_coordinate
    else:
        resolved_coordinate = coordinate

    adatom_point_defect_config = PointDefectSite(
        crystal=added_slab,
        element=AtomSchema(chemical_element=element),
        coordinate=resolved_coordinate,
    )
    isolated_defect = PointDefectSiteBuilder().get_material(adatom_point_defect_config)

    vacuum = slab_analyzer.get_slab_vacuum_configuration()

    configuration = AdatomDefectConfiguration(
        stack_components=[slab_without_vacuum, isolated_defect, vacuum],
    )

    builder = AdatomDefectBuilder()
    return builder.get_material(configuration)
