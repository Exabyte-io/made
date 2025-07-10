from typing import Optional, List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category_components.entities.auxiliary.zero_dimensional.point_defect_site import (
    AtomSchema,
)
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum
from sympy import ceiling

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.crystal_site import CrystalSiteAnalyzer, AdatomCrystalSiteAnalyzer
from mat3ra.made.tools.analyze.other import get_atomic_coordinates_extremum
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.build.defect.enums import AdatomPlacementMethodEnum
from mat3ra.made.tools.build.defect.point.builders import PointDefectSiteBuilder
from mat3ra.made.tools.build.defect.point.configuration import PointDefectSite
from mat3ra.made.tools.build.defect.slab.builders import AdatomDefectBuilder
from mat3ra.made.tools.build.defect.slab.configuration import AdatomDefectConfiguration
from ...defect.slab.builders import SlabStackBuilder
from ...defect.slab.configuration import SlabStackConfiguration
from ...slab.helpers import create_slab
from ....modify import filter_by_box


def create_slab_stack(slab: Material, added_component: Material) -> Material:
    """
    Create a slab stack by stacking a slab, a slab component, and a vacuum layer.

    Args:
        slab: The original slab material.
        added_component: The material to be stacked on top of the slab.
    Returns:
        Material: The new slab stack material.
    """
    analyzer = SlabMaterialAnalyzer(material=slab)

    slab_without_vacuum = analyzer.get_slab_configuration_with_no_vacuum()

    vacuum_config = analyzer.get_slab_vacuum_configuration()

    slab_stack_config = SlabStackConfiguration(
        stack_components=[slab_without_vacuum, added_component, vacuum_config], direction=AxisEnum.z
    )

    slab_stack_builder = SlabStackBuilder()
    return slab_stack_builder.get_material(slab_stack_config)


def recreate_slab_with_fractional_layers(slab: Material, number_of_layers: float) -> Material:
    """
    Create a slab with a specified number of fractional layers.

    Args:
        slab: The original slab material.
        number_of_layers: The total number of layers in the new slab.
    Returns:
        Material: The new slab material with the specified number of layers and vacuum if needed.
    """
    analyzer = SlabMaterialAnalyzer(material=slab)
    slab_without_vacuum = analyzer.get_slab_configuration_with_no_vacuum()
    # vacuum_config = analyzer.get_slab_vacuum_configuration()

    ceiling_number_of_layers = int(ceiling(number_of_layers))
    slab_with_int_layers_without_vacuum = create_slab(
        crystal=slab_without_vacuum.atomic_layers.crystal,
        miller_indices=slab_without_vacuum.atomic_layers.miller_indices,
        termination=slab_without_vacuum.atomic_layers.termination_top,
        number_of_layers=ceiling_number_of_layers,
        vacuum=0,
    )

    max_z_crystal_coordinate = number_of_layers / ceiling_number_of_layers
    return filter_by_box(
        slab_with_int_layers_without_vacuum,
        max_coordinate=[1, 1, max_z_crystal_coordinate],
        reset_ids=True,
    )


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
