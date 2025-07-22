from typing import List, Optional

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
        material=slab, coordinate_2d=position_on_surface, distance_z=distance_z, element=element
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
        material=slab, coordinate_2d=position_on_surface, distance_z=distance_z, element=element
    )

    slab_in_stack = adatom_analyzer.slab_configuration_with_no_vacuum
    isolated_defect = adatom_analyzer.added_component
    vacuum_configuration = adatom_analyzer.get_slab_vacuum_configuration()

    configuration = AdatomDefectConfiguration(
        stack_components=[slab_in_stack, isolated_defect, vacuum_configuration],
    )

    builder = AdatomDefectBuilder()
    return builder.get_material(configuration)


def create_multiple_adatom_defects(
    slab: MaterialWithBuildMetadata,
    adatom_dicts: List[dict],
    placement_method: AdatomPlacementMethodEnum = AdatomPlacementMethodEnum.EXACT_COORDINATE,
) -> Material:
    """
    Create multiple adatom defects from a list of dictionaries.

    Args:
        slab: The slab material.
        adatom_dicts: List of adatom dictionaries with keys:
            - element: str (chemical element for the adatom)
            - coordinate: List[float] (position on surface [x, y])
            - distance_z: float (distance above surface in Angstroms)
        placement_method: Method to place all adatoms (common for all defects).

    Returns:
        Material: The slab with all adatom defects applied.
    """
    all_adatom_configs = []

    for adatom_dict in adatom_dicts:
        element = adatom_dict["element"]
        coordinate = adatom_dict["coordinate"]
        distance_z = adatom_dict.get("distance_z", 1.0)

        if placement_method == AdatomPlacementMethodEnum.EXACT_COORDINATE:
            adatom_analyzer = AdatomMaterialAnalyzer(
                material=slab, coordinate_2d=coordinate, distance_z=distance_z, element=element
            )
        elif placement_method == AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE:
            adatom_analyzer = AdatomCrystalSiteMaterialAnalyzer(
                material=slab, coordinate_2d=coordinate, distance_z=distance_z, element=element
            )
        else:
            raise ValueError(f"Unsupported adatom placement method: {placement_method}")

        isolated_defect = adatom_analyzer.added_component
        all_adatom_configs.append(isolated_defect)

    # Create final configuration with all adatoms
    if placement_method == AdatomPlacementMethodEnum.EXACT_COORDINATE:
        slab_in_stack = adjust_material_cell_to_set_gap_along_direction(slab, 0)
        vacuum_configuration = AdatomMaterialAnalyzer(
            material=slab,
            coordinate_2d=adatom_dicts[0]["coordinate"],
            distance_z=adatom_dicts[0].get("distance_z", 1.0),
            element=adatom_dicts[0]["element"],
        ).get_slab_vacuum_configuration()
    else:
        analyzer = AdatomCrystalSiteMaterialAnalyzer(
            material=slab,
            coordinate_2d=adatom_dicts[0]["coordinate"],
            distance_z=adatom_dicts[0].get("distance_z", 1.0),
            element=adatom_dicts[0]["element"],
        )
        slab_in_stack = analyzer.slab_configuration_with_no_vacuum
        vacuum_configuration = analyzer.get_slab_vacuum_configuration()

    # Stack all components together
    stack_components = [slab_in_stack] + all_adatom_configs + [vacuum_configuration]
    configuration = AdatomDefectConfiguration(stack_components=stack_components)

    builder = AdatomDefectBuilder()
    return builder.get_material(configuration)
