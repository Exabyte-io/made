import numpy as np
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.made.material import Material

from ..core.unary import strain


def strain_to_match_lattice(material_to_strain: Material, target_material: Material) -> Material:
    """
    Strains a material to match the lattice of a target material.

    Args:
        material_to_strain: The material to be strained.
        target_material: The material with the target lattice.

    Returns:
        A new material instance with the strained lattice.
    """
    target_lattice_vectors = np.array(target_material.lattice.vector_arrays)
    source_lattice_vectors = np.array(material_to_strain.lattice.vector_arrays)

    # Prevent division by zero or singular matrix issues
    if np.linalg.det(source_lattice_vectors) == 0:
        raise ValueError("Lattice of material to strain is singular and cannot be inverted.")

    strain_matrix_np = np.linalg.inv(source_lattice_vectors) @ target_lattice_vectors
    strain_matrix = Matrix3x3Schema(root=strain_matrix_np.tolist())
    strained_material = strain(material_to_strain, strain_matrix)
    return strained_material
