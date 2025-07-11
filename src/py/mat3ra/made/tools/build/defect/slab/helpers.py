from typing import Union

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from sympy import ceiling

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.utils.coordinate import CoordinateCondition
from .builders import SlabStackBuilder
from .configuration import SlabStackConfiguration
from ...slab.helpers import create_slab
from ....modify import filter_by_box, filter_by_condition_on_coordinates


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


def create_island_defect(
    slab: Material,
    condition: CoordinateCondition,
    use_cartesian_coordinates: bool = True,
    number_of_added_layers: Union[int, float] = 1,
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

    material_with_additional_layers = recreate_slab_with_fractional_layers(slab, number_of_added_layers)

    isolated_island = filter_by_condition_on_coordinates(
        material=material_with_additional_layers,
        condition=condition.condition,
        use_cartesian_coordinates=use_cartesian_coordinates,
    )

    result = create_slab_stack(slab, isolated_island)
    return result
