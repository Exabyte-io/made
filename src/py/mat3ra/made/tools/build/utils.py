from typing import List, Optional, Union

from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.modify import translate_by_vector
from . import MaterialWithBuildMetadata
from .supercell import create_supercell
from ..modify import filter_by_box
from ..operations.core.binary import merge


def double_and_filter_material(
    material: Union[Material, MaterialWithBuildMetadata], start: List[float], end: List[float]
) -> Material:
    """
    Double the material and filter it by a box defined by the start and end coordinates.
    Args:
        material (Material): The material to double and filter.
        start (List[float]): The start coordinates of the box.
        end (List[float]): The end coordinates of the box.
    Returns:
        Material: The filtered material.
    """
    material_doubled = create_supercell(material, scaling_factor=[2, 1, 1])
    return filter_by_box(material_doubled, start, end)


def expand_lattice_vectors(
    material: Union[Material, MaterialWithBuildMetadata], gap: float, direction: int = 0
) -> Material:
    """
    Expand the lattice vectors of the material in the specified direction by the given gap.

    Args:
        material (Material): The material whose lattice vectors are to be expanded.
        gap (float): The gap by which to expand the lattice vector.
        direction (int): The index of the lattice vector to expand (0, 1, or 2).
    """
    new_lattice_vectors = material.lattice.vector_arrays
    new_lattice_vectors[direction][direction] += gap
    material.set_lattice_vectors(
        lattice_vector1=new_lattice_vectors[0],
        lattice_vector2=new_lattice_vectors[1],
        lattice_vector3=new_lattice_vectors[2],
    )
    return material


def stack_two_materials_xy(
    phase_1_material: Union[Material, MaterialWithBuildMetadata],
    phase_2_material: Union[Material, MaterialWithBuildMetadata],
    gap: float,
    edge_inclusion_tolerance: Optional[float] = 1.0,
    distance_tolerance: float = 1.0,
) -> Material:
    """
    Stack two materials laterally with translation along x-axis with a gap between them.

    Works correctly only for materials with the same lattice vectors (commensurate lattices).
     Args:
        phase_1_material (Material): The first material.
        phase_2_material (Material): The second material.
        gap (float): The gap between the two materials, in angstroms.
        edge_inclusion_tolerance (float): The tolerance to include atoms on the edge of the phase, in angstroms.
        distance_tolerance (float): The distance tolerance to remove atoms that are too close, in angstroms.

    Returns:
        Material: The merged material.
    """
    edge_inclusion_tolerance_crystal = abs(
        phase_1_material.basis.cell.convert_point_to_crystal([edge_inclusion_tolerance, 0, 0])[0]
    )

    phase_1_material = double_and_filter_material(
        phase_1_material, [0 - edge_inclusion_tolerance_crystal, 0, 0], [0.5 + edge_inclusion_tolerance_crystal, 1, 1]
    )

    phase_2_material = double_and_filter_material(
        phase_2_material, [0.5 - edge_inclusion_tolerance_crystal, 0, 0], [1 + edge_inclusion_tolerance_crystal, 1, 1]
    )

    phase_1_material = expand_lattice_vectors(phase_1_material, gap)
    phase_2_material = expand_lattice_vectors(phase_2_material, gap)

    phase_2_material = translate_by_vector(phase_2_material, [gap / 2, 0, 0], use_cartesian_coordinates=True)
    interface = merge(
        [phase_1_material, phase_2_material],
        merge_method=MergeMethodsEnum.ADD,
        distance_tolerance=distance_tolerance,
        merge_dangerously=True,
    )
    return interface
