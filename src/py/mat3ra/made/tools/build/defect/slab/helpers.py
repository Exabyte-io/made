from typing import Optional, List

from mat3ra.esse.models.materials_category_components.entities.auxiliary.zero_dimensional.point_defect_site import (
    AtomSchema,
)
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.crystal_site import CrystalSiteAnalyzer, AdatomCrystalSiteAnalyzer
from mat3ra.made.tools.analyze.other import get_atomic_coordinates_extremum
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.build.defect.configuration import SlabDefectConfigurationLegacy
from mat3ra.made.tools.build.defect.enums import AdatomPlacementMethodEnum
from mat3ra.made.tools.build.defect.point.builders import PointDefectSiteBuilder
from mat3ra.made.tools.build.defect.point.configuration import PointDefectSite
from mat3ra.made.tools.build.defect.slab.builders import SlabDefectBuilder, AdatomDefectBuilder, IslandDefectBuilder
from mat3ra.made.tools.build.defect.slab.configuration import SlabDefectConfiguration, AdatomDefectConfiguration, IslandDefectConfiguration, VoidSite
from mat3ra.made.tools.build.slab.builders import SlabBuilder
from mat3ra.made.tools.utils.coordinate import CoordinateCondition


def create_slab_defect(
    slab: Material, isolated_defect: Material, additional_layers: int = 1, vacuum_thickness: float = 0.0
) -> Material:
    """
    Create a slab defect by merging an isolated defect with a slab material.
    Args:
        slab: The original slab material.
        isolated_defect: The isolated defect material.
        additional_layers: Number of additional layers to add to the slab for defect creation.
    Returns:
        Material: The new slab material with additional layers and vacuum if needed.
    """
    analyzer = SlabMaterialAnalyzer(material=slab)
    (
        slab_with_additional_layers_config,
        slab_with_original_layers_config,
    ) = analyzer.get_slab_with_additional_layers_configuration_holder(
        additional_layers=additional_layers, vacuum_thickness=vacuum_thickness
    )

    # slab_with_additional_layers = SlabWithAdditionalLayersBuilder().get_material(slab_with_additional_layers_config)
    slab_with_original_layers = SlabBuilder().get_material(slab_with_original_layers_config)

    # defect is build with slab_with_additional_layers
    builder = SlabDefectBuilder()
    configuration = SlabDefectConfiguration(
        merge_components=[slab_with_original_layers, isolated_defect],
        merge_method=SlabDefectConfigurationLegacy.MergeMethodsEnum.ADD,
    )
    slab_with_defect = builder.get_material(configuration)
    return slab_with_defect


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
        analyzer = AdatomCrystalSiteAnalyzer(material=slab, coordinate=coordinate)
        resolved_coordinate = analyzer.new_crystal_site_coordinate
        slab = analyzer.get_slab_with_adjusted_vacuum()
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


def create_island_defect(
    slab: Material,
    condition: CoordinateCondition,
    number_of_added_layers: int = 1,
) -> Material:
    """
    Create an island defect using the new IslandDefectConfiguration and IslandDefectBuilder.
    
    Args:
        slab: The slab material.
        condition: The coordinate condition that defines the island shape.
        number_of_added_layers: Number of additional layers to add to the slab.
        
    Returns:
        Material: The slab with island defect.
    """
    # Create a slab with additional layers
    analyzer = SlabMaterialAnalyzer(material=slab)
    slab_with_additional_layers_config = analyzer.get_slab_with_additional_layers_configuration_holder(
        additional_layers=number_of_added_layers
    ).slab_with_additional_layers
    
    slab_with_additional_layers = SlabBuilder().get_material(slab_with_additional_layers_config)
    
    # Create a void site with the condition
    void_site = VoidSite(
        crystal=slab_with_additional_layers,
        condition=condition,
    )
    
    # Create the island configuration
    configuration = IslandDefectConfiguration(
        merge_components=[slab, void_site],
        merge_method=MergeMethodsEnum.REPLACE,
    )
    
    # Build the island defect
    builder = IslandDefectBuilder()
    return builder.get_material(configuration)
