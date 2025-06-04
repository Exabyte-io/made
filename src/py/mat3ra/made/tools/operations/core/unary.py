from typing import List

from mat3ra.code.vector import Vector3D

from mat3ra.made.material import Material
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.build.utils import stack_two_materials


def translate(material: Material, vector: Vector3D) -> Material:
    # reuse translate_by_vector, remove ASE
    # Figure out convention for use_cartesian_coordinates
    return material


def supercell(material: Material, supercell_matrix) -> Material:
    supercell = create_supercell(material, supercell_matrix)
    return supercell


def orient_cell(material: Material, rotational_matrix) -> Material:
    # reuse ase.rotate(rotate_cell=True)
    return material


def edit_cell(material: Material, lattice_vectors=None, basis=None) -> Material:
    # set lattice vectors, keep basis in cartesian
    return material


def stack(materials: List[Material], direction: AxisEnum) -> Material:
    # use stack_two_materials
    return stack_two_materials(material_1=materials[0], material_2=materials[1], direction=direction)
