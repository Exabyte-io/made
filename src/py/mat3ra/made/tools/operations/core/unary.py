from typing import List

from mat3ra.code.vector import Vector3D
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.build.utils import stack_two_materials
from mat3ra.made.tools.modify import translate_by_vector, wrap_to_unit_cell


def translate(material: Material, vector: Vector3D) -> Material:
    # Figure out convention for use_cartesian_coordinates
    return translate_by_vector(material, vector)


def supercell(material: Material, supercell_matrix) -> Material:
    supercell = create_supercell(material, supercell_matrix)
    return supercell


def edit_cell(material: Material, lattice_vectors=None) -> Material:
    # set lattice vectors, keep basis in cartesian
    if lattice_vectors is not None:
        material.set_new_lattice_vectors(
            lattice_vector1=lattice_vectors[0], lattice_vector2=lattice_vectors[1], lattice_vector3=lattice_vectors[2]
        )
    wrapped_material = wrap_to_unit_cell(material)
    return wrapped_material


def stack(materials: List[Material], direction: AxisEnum) -> Material:
    # use stack_two_materials
    return stack_two_materials(material_1=materials[0], material_2=materials[1], direction=direction)
