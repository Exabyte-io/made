from typing import Optional, List

from .build_parameters import TerraceBuildParameters
from .builder import TerraceDefectBuilder
from .configuration import TerraceDefectConfiguration
from .....analyze.terrace import TerraceMaterialAnalyzer
from .....build_components import MaterialWithBuildMetadata
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.builder import SlabBuilder
from .....build_components.entities.reusable.two_dimensional.slab_stack.helpers import (
    recreate_slab_with_fractional_layers,
)
from .....entities.coordinate import PlaneCoordinateCondition
from .....modify import filter_by_condition_on_coordinates


def create_defect_terrace(
    slab: MaterialWithBuildMetadata,
    cut_direction: Optional[List[int]] = None,
    pivot_coordinate: Optional[List[float]] = None,
    number_of_added_layers: float = 1.0,
    use_cartesian_coordinates: bool = False,
    rotate_to_match_pbc: bool = True,
) -> MaterialWithBuildMetadata:
    """
    Create a terrace at the specified position on the surface of the material.

    Args:
        slab: The original slab material to which the terrace will be added.
        cut_direction: The direction of the cut in lattice directions, Miller indices.
        pivot_coordinate: The center position of the terrace.
        number_of_added_layers: The number of added layers to the slab which will form the terrace
        use_cartesian_coordinates: Whether to use Cartesian coordinates for the center position.
        rotate_to_match_pbc: Whether to rotate the material to match the periodic boundary conditions.
    Returns:
        The material with the terrace added.
    """
    if cut_direction is None:
        cut_direction = [0, 0, 1]
    if pivot_coordinate is None:
        pivot_coordinate = [0.5, 0.5, 0.5]

    terrace_analyzer = TerraceMaterialAnalyzer(
        material=slab, cut_direction=cut_direction, number_of_added_layers=number_of_added_layers
    )

    additional_slab_material = recreate_slab_with_fractional_layers(slab, number_of_added_layers)

    normalized_direction_vector = terrace_analyzer.cut_direction_vector

    condition = PlaneCoordinateCondition(
        plane_normal=normalized_direction_vector,
        plane_point_coordinate=pivot_coordinate,
    ).condition

    slab_in_stack_config = terrace_analyzer.slab_configuration_with_no_vacuum
    slab_build_parameters = terrace_analyzer.build_parameters

    slab_in_stack = SlabBuilder(build_parameters=slab_build_parameters).get_material(slab_in_stack_config)
    isolated_defect = filter_by_condition_on_coordinates(
        material=additional_slab_material,
        condition=condition,
        use_cartesian_coordinates=use_cartesian_coordinates,
    )
    vacuum = terrace_analyzer.get_slab_vacuum_configuration()

    terrace_configuration = TerraceDefectConfiguration(
        stack_components=[slab_in_stack, isolated_defect, vacuum], cut_direction=cut_direction
    )

    angle = terrace_analyzer.angle
    axis = terrace_analyzer.rotation_axis
    new_lattice_vectors = terrace_analyzer.new_lattice_vectors
    parameters = TerraceBuildParameters(
        rotate_to_match_pbc=rotate_to_match_pbc, angle=angle, axis=axis, new_lattice_vectors=new_lattice_vectors
    )
    terrace = TerraceDefectBuilder(build_parameters=parameters).get_material(terrace_configuration)

    return terrace
