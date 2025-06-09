import numpy as np

from mat3ra.made.material import Material
from ...operations.core.unary import edit_cell


def get_orthogonal_c_slab(material: Material) -> Material:
    """
    Make the c-vector orthogonal to the ab plane after slab construction is complete.
    This should be applied after vacuum has been added to the slab.
    """
    current_vectors = material.lattice.vector_arrays
    a = np.array(current_vectors[0])
    b = np.array(current_vectors[1])
    c_old = np.array(current_vectors[2])

    normal = np.cross(a, b)
    normal_unit = np.linalg.norm(normal)
    n_hat = normal / normal_unit

    height = float(np.dot(c_old, n_hat))

    new_orthogonal_vector_c = (n_hat * height).tolist()

    new_vectors = [
        current_vectors[0],
        current_vectors[1],
        new_orthogonal_vector_c,
    ]

    return edit_cell(material, new_vectors)
