from types import SimpleNamespace
from typing import List, Optional

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.adatom import AdatomCrystalSiteMaterialAnalyzer, AdatomMaterialAnalyzer
from mat3ra.made.tools.operations.core.binary import merge
from .builders import AdatomDefectBuilder
from .configuration import (
    AdatomDefectConfiguration,
)
from ... import MaterialWithBuildMetadata
from ...defect.enums import AdatomPlacementMethodEnum


def get_adatom_defect_analyzer_cls(
    placement_method: str = AdatomPlacementMethodEnum.EXACT_COORDINATE.value,
):
    if placement_method == AdatomPlacementMethodEnum.EXACT_COORDINATE.value:
        return AdatomMaterialAnalyzer
    elif placement_method == AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE.value:
        return AdatomCrystalSiteMaterialAnalyzer
    else:
        raise ValueError(f"Unsupported placement method: {placement_method}")


def create_adatom_defect(
    slab: MaterialWithBuildMetadata,
    position_on_surface: List[float],
    distance_z: float = 1.0,
    placement_method: str = AdatomPlacementMethodEnum.EXACT_COORDINATE.value,
    element: Optional[str] = None,
    use_cartesian_coordinates: bool = False,
) -> Material:
    """
    Create an adatom defect based on the specified placement method.

    Args:
        slab: The slab material.
        position_on_surface: Position on the surface [x, y].
        distance_z: Distance above the surface in Angstroms.
        placement_method: Method to place the adatom.
        element: Chemical element for the adatom.
        use_cartesian_coordinates: Whether the position_on_surface is in Cartesian units.

    Returns:
        Material: The slab with adatom defect.
    """
    # Reuse the create_multiple_adatom_defects function below
    slab_with_adatom = create_multiple_adatom_defects(
        slab,
        defect_dicts=[
            {
                "element": element or "Si",  # Default to Silicon if no element provided
                "coordinate": position_on_surface,
                "distance_z": distance_z,
                "use_cartesian_coordinates": use_cartesian_coordinates,
            }
        ],
        placement_method=placement_method,
    )
    return slab_with_adatom


def create_multiple_adatom_defects(
    slab: MaterialWithBuildMetadata, defect_dicts: List[dict], placement_method: str
) -> Material:
    """
    Create multiple adatom defects (at once) from a list of dictionaries.

    Args:
        slab: The slab material.
        defect_dicts: List of adatom dictionaries with keys:
            - element: str (chemical element for the adatom)
            - coordinate: List[float] (position on surface [x, y])
            - distance_z: float (distance above surface in Angstroms)
            - use_cartesian_coordinates: bool (optional, defaults to False)
                Whether the coordinate is in Cartesian units.
        placement_method: Method to place all adatoms (common for all defects).
            Valid values from AtomPlacementMethodEnum:
                - "exact_coordinate" (default): Places atom at exact coordinate
                - "new_crystal_site": Places atom at nearest crystal site
                - "equidistant": Places atom equidistant from nearest atoms

    Returns:
        Material: The slab with all adatom defects applied.
    """
    if placement_method not in [e.value for e in AdatomPlacementMethodEnum]:
        raise ValueError(f"Unsupported placement method: {placement_method}")

    all_adatom_configs = []
    analyzer_cls = get_adatom_defect_analyzer_cls(placement_method)

    last_analyzer = None

    for defect_dict in defect_dicts:
        defect_configuration = SimpleNamespace(**defect_dict)

        coordinate_2d = defect_configuration.coordinate
        use_cartesian = getattr(defect_configuration, "use_cartesian_coordinates", False)

        if use_cartesian:
            coordinate_3d = coordinate_2d + [0.0]
            coordinate_3d_crystal = slab.basis.cell.convert_point_to_crystal(coordinate_3d)
            coordinate_2d = coordinate_3d_crystal[:2]

        analyzer = analyzer_cls(
            material=slab,
            coordinate_2d=coordinate_2d,
            distance_z=defect_configuration.distance_z,
            placement_method=placement_method,
            element=defect_configuration.element,
        )
        last_analyzer = analyzer
        all_adatom_configs.append(analyzer.added_component)

    vacuum_configuration = last_analyzer.get_slab_vacuum_configuration()

    merged_adatoms_component = merge(all_adatom_configs)

    stack_components = [
        last_analyzer.slab_material_or_configuration_for_stacking,
        merged_adatoms_component,
        vacuum_configuration,
    ]
    configuration = AdatomDefectConfiguration(stack_components=stack_components)

    builder = AdatomDefectBuilder()
    return builder.get_material(configuration)
