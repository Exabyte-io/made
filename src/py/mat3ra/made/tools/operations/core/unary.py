from typing import List

from mat3ra.code.vector import Vector3D

from mat3ra.made.material import Material


def translate(material: Material, vector: Vector3D) -> Material:
    # reuse translate_by_vector, remove ASE
    # Figure out convention for use_cartesian_coordinates
    pass


def supercell(material: Material, supercell_matrix) -> Material:
    # reuse create_supercell, keep ASE
    pass


def orient_cell(material: Material, rotational_matrix) -> Material:
    # reuse ase.rotate(rotate_cell=True)
    pass


def edit_cell():
    # set lattice vectors, keep basis in cartesian
    pass


def stack(materials: List[Material], direction: str) -> Material:
    # use stack_two_materials
    pass
