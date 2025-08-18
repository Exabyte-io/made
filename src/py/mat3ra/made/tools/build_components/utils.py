from typing import List, Union

from mat3ra.made.material import Material

from ..modify import filter_by_box
from . import MaterialWithBuildMetadata
from .entities.reusable.three_dimensional.supercell.helpers import create_supercell


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
