from typing import List, Optional

from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import (
    AtomSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.adatom import AdatomCrystalSiteMaterialAnalyzer, AdatomMaterialAnalyzer
from mat3ra.made.utils import adjust_material_cell_to_set_gap_along_direction
from .builders import AdatomDefectBuilder
from .configuration import (
    AdatomDefectConfiguration,
)
from ... import MaterialWithBuildMetadata
from ...defect.enums import AdatomPlacementMethodEnum


def create_adatom_defect(
    slab: MaterialWithBuildMetadata,
    position_on_surface: List[float],
    distance_z: float = 1.0,
    placement_method: AdatomPlacementMethodEnum = AdatomPlacementMethodEnum.EXACT_COORDINATE,
    element: Optional[str] = None,
) -> Material:
    """
    Create an adatom defect based on the specified placement method.

    Args:
        slab: The slab material.
        position_on_surface: Position on the surface [x, y].
        distance_z: Distance above the surface in Angstroms.
        placement_method: Method to place the adatom.
        element: Chemical element for the adatom.

    Returns:
        Material: The slab with adatom defect.
    """

    if placement_method == AdatomPlacementMethodEnum.EXACT_COORDINATE:
        return create_adatom_defect_at_exact_coordinate(slab, position_on_surface, distance_z, element)
    elif placement_method == AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE:
        return create_adatom_defect_at_crystal_site(slab, position_on_surface, distance_z, element)
    else:
        raise ValueError(f"Unsupported adatom placement method: {placement_method}")


def create_adatom_defect_at_exact_coordinate(
    slab: MaterialWithBuildMetadata,
    position_on_surface: List[float],
    distance_z: float = 1.0,
    element: Optional[str] = None,
) -> Material:
    """
    Create an adatom defect using the new AdatomDefectConfiguration and AdatomDefectBuilder.

    Args:
        slab: The slab material.
        position_on_surface: Position on the surface [x, y].
        element: Chemical element for the adatom.
        distance_z: Distance above the surface in Angstroms.

    Returns:
        Material: The slab with adatom defect.
    """

    adatom_analyzer = AdatomMaterialAnalyzer(
        material=slab,
        coordinate_2d=position_on_surface,
        distance_z=distance_z,
        element=AtomSchema(chemical_element=element),
    )

    slab_in_stack = adjust_material_cell_to_set_gap_along_direction(slab, 0)
    isolated_defect = adatom_analyzer.added_component
    vacuum_configuration = adatom_analyzer.get_slab_vacuum_configuration()

    configuration = AdatomDefectConfiguration(
        stack_components=[slab_in_stack, isolated_defect, vacuum_configuration],
    )

    builder = AdatomDefectBuilder()
    return builder.get_material(configuration)


def create_adatom_defect_at_crystal_site(
    slab: MaterialWithBuildMetadata,
    position_on_surface: List[float],
    distance_z: float = 1.0,
    element: Optional[str] = None,
) -> Material:
    """
    Create an adatom defect at a crystal site.

    Args:
        slab: The slab material.
        position_on_surface: Position on the surface [x, y].
        element: Chemical element for the adatom.
        distance_z: Distance above the surface in Angstroms.

    Returns:
        Material: The slab with adatom defect.
    """

    adatom_analyzer = AdatomCrystalSiteMaterialAnalyzer(
        material=slab,
        coordinate_2d=position_on_surface,
        distance_z=distance_z,
        element=AtomSchema(chemical_element=element),
    )

    slab_in_stack = adatom_analyzer.slab_configuration_with_no_vacuum
    isolated_defect = adatom_analyzer.added_component
    vacuum_configuration = adatom_analyzer.get_slab_vacuum_configuration()

    configuration = AdatomDefectConfiguration(
        stack_components=[slab_in_stack, isolated_defect, vacuum_configuration],
    )

    builder = AdatomDefectBuilder()
    return builder.get_material(configuration)
