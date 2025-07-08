from mat3ra.esse.models.materials_category_components.entities.auxiliary.zero_dimensional.point_defect_site import (
    AtomSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.crystal_site import CrystalSiteAnalyzer
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.build.defect.configuration import SlabDefectConfigurationLegacy
from mat3ra.made.tools.build.defect.enums import AdatomPlacementMethodEnum
from mat3ra.made.tools.build.defect.point.builders import PointDefectBuilder
from mat3ra.made.tools.build.defect.point.configuration import PointDefectSite
from mat3ra.made.tools.build.defect.slab.builders import SlabDefectBuilder
from mat3ra.made.tools.build.defect.slab.configuration import SlabDefectConfiguration
from mat3ra.made.tools.build.slab.builders import SlabBuilder


def create_slab_defect(slab: Material, isolated_defect: Material, additional_layers: int = 1) -> Material:
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
    ) = analyzer.get_slab_with_additional_layers_configurations(
        additional_layers=additional_layers, vacuum_thickness=5.0
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
    position_on_surface: list[float],
    placement_method: AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE,
    element: str,
    added_vacuum: float = 5.0,
) -> Material:
    max_z = slab.basis.coordinates.get_extremum_value_along_axis("max", "z")
    coordinate = [*position_on_surface, max_z]
    if placement_method == AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE:
        analyzer = SlabMaterialAnalyzer(material=slab)
        (
            slab_with_additional_layers_config,
            slab_with_original_layers_config,
        ) = analyzer.get_slab_with_additional_layers_configurations(additional_layers=1, vacuum_thickness=added_vacuum)
        slab = SlabBuilder().get_material(slab_with_original_layers_config)
        slab_with_added_layer = SlabBuilder().get_material(slab_with_additional_layers_config)
        crystal_site_analyzer = CrystalSiteAnalyzer(material=slab_with_added_layer, coordinate=coordinate)
        resolved_coordinate = crystal_site_analyzer.closest_site_coordinate

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

    isolated_defect = PointDefectBuilder().get_material(adatom_point_defect_config)

    configuration = SlabDefectConfiguration.from_materials(
        slab=slab,
        isolated_defect=isolated_defect,
        merge_method=SlabDefectConfiguration.MergeMethodsEnum.ADD,
    )

    builder = SlabDefectBuilder()
    slab_with_adatom_defect = builder.get_material(configuration)
    return slab_with_adatom_defect
