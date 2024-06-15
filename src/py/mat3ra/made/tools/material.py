import numpy as np
from mat3ra.code.constants import AtomicCoordinateUnits


def scale_one_lattice_vector(material, key="a", factor=1.0):
    """
    Scales one lattice vector for the given material.
    :param material: The material acted upon.
    :param key: Lattice vector key.
    :param factor: Float scaling factor.
    """
    material.to_cartesian()

    lattice = getattr(material, "lattice")
    setattr(lattice, key, getattr(lattice, key) * factor)

    setattr(material, "lattice", lattice)

    material.to_crystal()


def get_basis_config_translated_to_center(material):
    """
    Updates the basis of a material by translating the coordinates
    so that the center of the material and lattice are aligned.
    :param material: The material to update.
    """
    original_units = material.basis.units
    material.to_cartesian()
    updated_basis = material.basis
    center_of_coordinates = updated_basis.center_of_coordinates_point()
    center_of_lattice = 0.5 * np.sum(material.lattice.vector_arrays, axis=0)
    translation_vector = center_of_lattice - center_of_coordinates
    updated_basis.translate_by_vector(translation_vector)
    material.set_basis(updated_basis.to_json())
    if original_units != AtomicCoordinateUnits.cartesian:
        material.to_crystal()
