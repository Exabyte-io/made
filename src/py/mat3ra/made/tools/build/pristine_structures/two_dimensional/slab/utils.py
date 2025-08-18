from typing import Union, Optional, List

import numpy as np

from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from .....build_components import MaterialWithBuildMetadata
from .....build_components.entities.auxiliary.two_dimensional.termination import Termination
from .....modify import wrap_to_unit_cell
from .....operations.core.unary import edit_cell


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
