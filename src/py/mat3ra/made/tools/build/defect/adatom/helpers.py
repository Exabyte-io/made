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


def create_adatom_defect(
    slab: Material,
    position_on_surface: List[float],
    distance_z: float = 1.0,
    placement_method: AdatomPlacementMethodEnum = AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE,
    element: Optional[str] = None,
    added_vacuum: float = 5.0,
) -> Material:
    """
    Create an adatom defect using the new AdatomDefectConfiguration and AdatomDefectBuilder.

    Args:
        slab: The slab material.
        position_on_surface: Position on the surface [x, y].
        element: Chemical element for the adatom.
        distance_z: Distance above the surface in Angstroms.
        placement_method: Method for placing the adatom.
        added_vacuum: Additional vacuum thickness.

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
        crystal_site_analyzer = CrystalSiteAnalyzer(material=added_slab, coordinate=coordinate)
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


#
#
# def create_crystal_site_adatom_defect(
#     slab: Material,
#     position_on_surface: List[float],
#     distance_z: float = 1.0,
#     element: Optional[str] = None,
# ) -> Material:
#     """
#     Create an adatom defect using crystal site placement with the new SlabStackBuilder approach.
#
#     This function uses the new approach where:
#     1. Original slab (without vacuum)
#     2. Added component (single atom in lattice from recreate_slab_with_fractional_layers)
#     3. Vacuum layer
#
#     Args:
#         slab: The slab material.
#         position_on_surface: Position on the surface [x, y].
#         distance_z: Distance above the surface in Angstroms.
#         element: Chemical element for the adatom.
#
#     Returns:
#         Material: The slab with adatom defect using crystal site approach.
#     """
#     distance_z_crystal = slab.basis.cell.convert_point_to_crystal([0, 0, distance_z])[2]
#     coordinate = [*position_on_surface, distance_z_crystal]
#     analyzer = AdatomCrystalSiteAnalyzer(material=slab, coordinate=coordinate)
#     coordinate = analyzer.new_crystal_site_coordinate
#
#     configuration = AdatomDefectConfiguration.from_parameters(slab=slab, coordinate=coordinate, element=element)
#
#     builder = AdatomDefectBuilder()
#     return builder.get_material(configuration)
