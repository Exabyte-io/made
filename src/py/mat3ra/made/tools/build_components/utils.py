from typing import List, Optional, Union

import numpy as np
from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material

from ..modify import filter_by_box, wrap_to_unit_cell
from ..operations.core.unary import edit_cell
from . import MaterialWithBuildMetadata
from .entities.auxiliary.two_dimensional.termination import Termination
from .entities.reusable.three_dimensional.supercell.helpers import create_supercell


def select_slab_termination(terminations: List[Termination], formula: Optional[str] = None) -> Termination:
    if not terminations:
        raise ValueError("No terminations available.")
    if formula is None:
        return terminations[0]
    for termination in terminations:
        if termination.formula == formula:
            return termination
    raise ValueError(f"Termination with formula {formula} not found in available terminations: {terminations}")


def get_orthogonal_c_slab(material: Union[Material, MaterialWithBuildMetadata]) -> Material:
    """
    Make the c-vector orthogonal to the ab plane and update the basis.

    This function calculates a new c-vector that is orthogonal to the a and b vectors
    of the lattice. It then computes the transformation matrix between the old and new
    lattice vectors and applies this transformation to the atomic coordinates.

    A new material is returned with an updated lattice and basis, where the new
    lattice is defined by its parameters (a, b, c, alpha, beta, gamma) to avoid
    storing raw vectors, preserving a standard representation.

    Args:
        material (Material): The input material object.

    Returns:
        Material: A new material object with an orthogonalized c-vector and
                  updated basis.
    """
    new_material = material.clone()
    current_vectors = np.array(new_material.lattice.vector_arrays)
    a_vec, b_vec, c_old_vec = current_vectors

    normal = np.cross(a_vec, b_vec)
    n_hat = normal / np.linalg.norm(normal)
    height = float(np.dot(c_old_vec, n_hat))
    c_new_vec = n_hat * height

    new_vectors = np.array([a_vec, b_vec, c_new_vec])
    transform_matrix = np.dot(current_vectors, np.linalg.inv(new_vectors))

    new_basis = new_material.basis.clone()
    new_basis.transform_by_matrix(transform_matrix)
    new_lattice_from_vectors = Lattice.from_vectors_array(new_vectors.tolist())
    new_material = edit_cell(new_material, new_lattice_from_vectors.vector_arrays)
    new_material.basis = new_basis
    new_material = wrap_to_unit_cell(new_material)
    return new_material


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
