from mat3ra.code.vector import Vector3D

from mat3ra.made.material import Material
from mat3ra.made.tools.modify import translate_by_vector, wrap_to_unit_cell
from ...convert import from_ase, to_ase
from ...third_party import ase_make_supercell
from ...utils import decorator_convert_supercell_matrix_2x2_to_3x3


def translate(material: Material, vector: Vector3D) -> Material:
    # Figure out convention for use_cartesian_coordinates
    return translate_by_vector(material, vector)


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
