from typing import List, Optional

from mat3ra.esse.models.materials_category.defective_structures.two_dimensional.adatom.configuration import AtomSchema
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.tools.build.defect.adatom.builders import AdatomDefectBuilder, CrystalSiteAdatomSlabDefectBuilder
from mat3ra.made.tools.build.defect.adatom.configuration import (
    AdatomDefectConfiguration,
    CrystalSiteAdatomConfiguration,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.crystal_site import CrystalSiteAnalyzer
from mat3ra.made.tools.analyze.other import get_atomic_coordinates_extremum
from mat3ra.made.tools.build.defect.enums import AdatomPlacementMethodEnum
from mat3ra.made.tools.build.defect.point.builders import PointDefectSiteBuilder
from mat3ra.made.tools.build.defect.point.configuration import PointDefectSite


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
    max_z = get_atomic_coordinates_extremum(slab, "max", "z")
    distance_z_crystal = slab.basis.cell.convert_point_to_crystal([0, 0, distance_z])[2]
    coordinate = [*position_on_surface, max_z + distance_z_crystal]

    if placement_method == AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE:
        return create_crystal_site_adatom_defect(
            slab=slab,
            position_on_surface=position_on_surface,
            distance_z=distance_z,
            element=element,
        )
    elif placement_method == AdatomPlacementMethodEnum.EQUIDISTANT:
        crystal_site_analyzer = CrystalSiteAnalyzer(material=slab, coordinate=coordinate)
        resolved_coordinate = crystal_site_analyzer.equidistant_coordinate
    else:
        resolved_coordinate = coordinate

    adatom_point_defect_config = PointDefectSite(
        crystal=slab,
        element=AtomSchema(chemical_element=element),
        coordinate=resolved_coordinate,
    )
    isolated_defect = PointDefectSiteBuilder().get_material(adatom_point_defect_config)

    configuration = AdatomDefectConfiguration(
        merge_components=[slab, isolated_defect],
        merge_method=MergeMethodsEnum.ADD,
    )

    builder = AdatomDefectBuilder()
    return builder.get_material(configuration)


def create_crystal_site_adatom_defect(
    slab: Material,
    position_on_surface: List[float],
    distance_z: float = 1.0,
    element: Optional[str] = None,
) -> Material:
    """
    Create an adatom defect using crystal site placement with the new SlabStackBuilder approach.

    This function uses the new approach where:
    1. Original slab (without vacuum)
    2. Added component (single atom in lattice from recreate_slab_with_fractional_layers)
    3. Vacuum layer

    Args:
        slab: The slab material.
        position_on_surface: Position on the surface [x, y].
        distance_z: Distance above the surface in Angstroms.
        element: Chemical element for the adatom.

    Returns:
        Material: The slab with adatom defect using crystal site approach.
    """
    distance_z_crystal = slab.basis.cell.convert_point_to_crystal([0, 0, distance_z])[2]
    coordinate = [*position_on_surface, distance_z_crystal]

    configuration = CrystalSiteAdatomConfiguration.from_parameters(slab=slab, coordinate=coordinate, element=element)

    builder = CrystalSiteAdatomSlabDefectBuilder()
    return builder.get_material(configuration)
