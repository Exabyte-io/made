from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from sympy import ceiling

from mat3ra.made.material import Material
from .builders import SlabStackBuilder
from .configuration import SlabStackConfiguration
from ...slab.helpers import create_slab
from ....analyze.slab import SlabMaterialAnalyzer
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

    slab_without_vacuum = analyzer.slab_configuration_with_no_vacuum

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
    slab_without_vacuum = analyzer.slab_configuration_with_no_vacuum
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
