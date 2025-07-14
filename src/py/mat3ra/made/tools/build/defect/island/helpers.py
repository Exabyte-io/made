from typing import Union

from mat3ra.made.material import Material
from ..slab.helpers import recreate_slab_with_fractional_layers, create_slab_stack
from ....modify import filter_by_condition_on_coordinates
from ....utils.coordinate import CoordinateCondition


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
