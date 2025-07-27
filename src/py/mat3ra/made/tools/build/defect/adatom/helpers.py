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

# Define mapping of placement methods to their corresponding analyzers
ADATOM_PLACEMENT_MAPPING = {
    AdatomPlacementMethodEnum.EXACT_COORDINATE: {
        "analyzer_class": AdatomMaterialAnalyzer,
        "needs_gap_adjustment": True,
    },
    AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE: {
        "analyzer_class": AdatomCrystalSiteMaterialAnalyzer,
        "needs_gap_adjustment": False,
    },
}


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
    # Reuse the create_multiple_adatom_defects function below
    slab_with_adatom = create_multiple_adatom_defects(
        slab,
        adatom_dicts=[
            {
                "element": element or "Si",  # Default to Silicon if no element provided
                "coordinate": position_on_surface,
                "distance_z": distance_z,
            }
        ],
        placement_method=placement_method,
    )
    return slab_with_adatom


def create_multiple_adatom_defects(
    slab: MaterialWithBuildMetadata,
    adatom_dicts: List[dict],
    placement_method: AdatomPlacementMethodEnum = AdatomPlacementMethodEnum.EXACT_COORDINATE,
) -> Material:
    """
    Create multiple adatom defects (at once) from a list of dictionaries.

    Args:
        slab: The slab material.
        adatom_dicts: List of adatom dictionaries with keys:
            - element: str (chemical element for the adatom)
            - coordinate: List[float] (position on surface [x, y])
            - distance_z: float (distance above surface in Angstroms)
        placement_method: Method to place all adatoms (common for all defects).
            Valid values from AtomPlacementMethodEnum:
                - "exact_coordinate" (default): Places atom at exact coordinate
                - "new_crystal_site": Places atom at nearest crystal site
                - "equidistant": Places atom equidistant from nearest atoms

    Returns:
        Material: The slab with all adatom defects applied.
    """
    if placement_method not in ADATOM_PLACEMENT_MAPPING:
        raise ValueError(f"Unsupported adatom placement method: {placement_method}")

    placement_info = ADATOM_PLACEMENT_MAPPING[placement_method]
    all_adatom_configs = []

    # Create all adatom components
    for adatom_dict in adatom_dicts:
        analyzer = placement_info["analyzer_class"](
            material=slab,
            coordinate_2d=adatom_dict["coordinate"],
            distance_z=adatom_dict.get("distance_z", 1.0),
            element=adatom_dict["element"],
        )
        all_adatom_configs.append(analyzer.added_component)

    # Create final configuration
    analyzer = placement_info["analyzer_class"](
        material=slab,
        coordinate_2d=adatom_dicts[0]["coordinate"],
        distance_z=adatom_dicts[0].get("distance_z", 1.0),
        element=adatom_dicts[0]["element"],
    )

    slab_in_stack = (
        adjust_material_cell_to_set_gap_along_direction(slab, 0)
        if placement_info["needs_gap_adjustment"]
        else analyzer.slab_configuration_with_no_vacuum
    )
    vacuum_configuration = analyzer.get_slab_vacuum_configuration()

    # Stack all components together
    stack_components = [slab_in_stack] + all_adatom_configs + [vacuum_configuration]
    configuration = AdatomDefectConfiguration(stack_components=stack_components)

    builder = AdatomDefectBuilder()
    return builder.get_material(configuration)
