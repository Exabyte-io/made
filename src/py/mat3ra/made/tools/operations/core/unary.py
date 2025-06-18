import numpy as np
from mat3ra.code.vector import Vector3D
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.modify import translate_by_vector, wrap_to_unit_cell

from ...convert import from_ase, to_ase
from ...third_party import ase_make_supercell
from ...utils import decorator_convert_supercell_matrix_2x2_to_3x3


def translate(material: Material, vector: Vector3D) -> Material:
    # Figure out convention for use_cartesian_coordinates
    return translate_by_vector(material, vector, use_cartesian_coordinates=True)


@decorator_convert_supercell_matrix_2x2_to_3x3
def supercell(material: Material, supercell_matrix) -> Material:
    atoms = to_ase(material)

    supercell_atoms = ase_make_supercell(atoms, supercell_matrix)
    new_material = Material.create(from_ase(supercell_atoms))
    if material.metadata:
        new_material.metadata = material.metadata
    new_material.name = material.name
    return new_material


def edit_cell(material: Material, lattice_vectors=None) -> Material:
    if lattice_vectors is not None:
        material.set_new_lattice_vectors(
            lattice_vector1=lattice_vectors[0], lattice_vector2=lattice_vectors[1], lattice_vector3=lattice_vectors[2]
        )
    wrapped_material = wrap_to_unit_cell(material)
    return wrapped_material


def strain(material: Material, strain_matrix: Matrix3x3Schema) -> Material:
    # Applies a strain to the material by modifying its lattice vectors while keeping the basis in crystal coordinates.

    # The strain matrix is a 3x3 matrix that defines the strain to be applied.
    # It should be a small perturbation to the identity matrix.
    if not isinstance(strain_matrix, Matrix3x3Schema):
        raise ValueError("strain_matrix must be an instance of Matrix3x3Schema")

    strain_matrix_np = np.array([row.root for row in strain_matrix.root])
    lattice_vectors_np = np.array(material.lattice.vector_arrays)
    new_lattice_vectors_np = lattice_vectors_np @ strain_matrix_np
    new_lattice_vectors = new_lattice_vectors_np.tolist()

    new_material = material.clone()

    original_crystal_coords = new_material.basis.coordinates.values

    new_material.set_new_lattice_vectors(
        lattice_vector1=new_lattice_vectors[0],
        lattice_vector2=new_lattice_vectors[1],
        lattice_vector3=new_lattice_vectors[2],
    )

    new_material.basis.coordinates.values = original_crystal_coords

    return new_material


def mirror(material: Material, direction: AxisEnum = AxisEnum.z) -> Material:
    # Flips the material along the specified axis by applying a right-handed supercell transformation.
    supercell_matrix = {
        AxisEnum.x: [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],
        AxisEnum.y: [[0, 0, 1], [0, -1, 0], [1, 0, 0]],
        AxisEnum.z: [[0, 1, 0], [1, 0, 0], [0, 0, -1]],
    }.get(direction)

    mirrored_material = supercell(material, supercell_matrix)
    return mirrored_material
