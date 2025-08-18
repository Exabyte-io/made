from collections import Counter
from typing import List, Dict

from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import AtomSchema

from .analyzers import (
    SurfacePassivationMaterialAnalyzer,
    CoordinationBasedPassivationMaterialAnalyzer,
)
from .builder import PassivationBuilder
from .configuration import PassivationConfiguration
from .... import MaterialWithBuildMetadata
from .....analyze.material import MaterialWithCrystalSites
from .....build_components.entities.core.zero_dimensional.atom.configuration import AtomAtCoordinateConfiguration
from mat3ra.made.tools.build.processed_structures.two_dimensional.passivation.enums import SurfaceTypesEnum


def passivate_surface(
    material: MaterialWithBuildMetadata,
    passivant: str = "H",
    bond_length: float = 1.0,
    shadowing_radius: float = 2.5,
    depth: float = 5.0,
    surface: SurfaceTypesEnum = SurfaceTypesEnum.TOP,
) -> MaterialWithBuildMetadata:
    """
    Create a passivated material by adding passivant atoms to the surface of the provided material.

    Args:
        material (Material): The material to be passivated. Should have vacuum space above or below the surface.
        passivant (str): The chemical symbol of the passivating atom (default is 'H').
        bond_length (float): The bond length for the passivation (default is 1.0). Å
        shadowing_radius (float): The radius with which an atom shadows/covers underlying atom from being exposed. Å
        depth (float): The depth from the most top (or bottom) atom into the slab to search for exposed atoms. Å
        surface (SurfaceTypesEnum): The surface to passivate (TOP, BOTTOM, or BOTH).
    Returns:
        Material: The passivated material.
    """
    analyzer = SurfacePassivationMaterialAnalyzer(
        material=material,
        passivant=passivant,
        bond_length=bond_length,
        shadowing_radius=shadowing_radius,
        depth=depth,
        surface=surface,
    )

    passivant_coordinates = analyzer.passivant_coordinates
    passivant_configs = [
        AtomAtCoordinateConfiguration(
            crystal=material, element=AtomSchema(chemical_element=passivant), coordinate=coord
        )
        for coord in passivant_coordinates
    ]
    config = PassivationConfiguration(
        merge_components=[material, *passivant_configs], passivant=passivant, bond_length=bond_length
    )

    builder = PassivationBuilder()
    return builder.get_material(config)


def passivate_dangling_bonds(
    material: MaterialWithBuildMetadata,
    passivant: str = "H",
    bond_length: float = 1.0,
    coordination_threshold: int = 3,
    number_of_bonds_to_passivate: int = 1,
    symmetry_tolerance: float = 0.1,
    shadowing_radius: float = 2.5,
    depth: float = 5.0,
) -> MaterialWithBuildMetadata:
    """
    Create a passivated material by adding passivant atoms to the undercoordinated atoms of the provided material.

    Args:
        material (MaterialWithBuildMetadata): The material to be passivated.
        passivant (str): The chemical symbol of the passivating atom (default is 'H').
        bond_length (float): The bond length for the passivation (default is 1.0). Å
        coordination_threshold (int): The coordination number threshold for an atom to be considered undercoordinated.
        number_of_bonds_to_passivate (int): The maximum number of bonds to passivate for each undercoordinated atom.
        symmetry_tolerance (float): The tolerance for symmetry comparison of vectors for bonds.
        shadowing_radius (float): The radius with which an atom shadows/covers underlying atom from being exposed. Å
        depth (float): The depth from the most top (or bottom) atom into the slab to search for exposed atoms. Å

    Returns:
        MaterialWithBuildMetadata: The passivated material with dangling bonds filled.
    """
    analyzer = CoordinationBasedPassivationMaterialAnalyzer(
        material=material,
        passivant=passivant,
        bond_length=bond_length,
        coordination_threshold=coordination_threshold,
        number_of_bonds_to_passivate=number_of_bonds_to_passivate,
        symmetry_tolerance=symmetry_tolerance,
        shadowing_radius=shadowing_radius,
        depth=depth,
    )

    passivant_coordinates = analyzer.passivant_coordinates
    passivant_configs = [
        AtomAtCoordinateConfiguration(
            crystal=material, element=AtomSchema(chemical_element=passivant), coordinate=coord
        )
        for coord in passivant_coordinates
    ]
    config = PassivationConfiguration(
        merge_components=[material, *passivant_configs], passivant=passivant, bond_length=bond_length
    )

    builder = PassivationBuilder()
    return builder.get_material(config)


def get_unique_coordination_numbers(
    material: MaterialWithBuildMetadata,
    cutoff: float = 3.0,
) -> List[int]:
    """
    Get the unique coordination numbers for the provided passivation configuration and cutoff radius.

    Args:
        material (MaterialWithBuildMetadata): The material object.
        cutoff (float): The cutoff radius for defining neighbors.
    Returns:
        set: The unique coordination numbers.
    """
    material_with_crystal_sites = MaterialWithCrystalSites.from_material(material)
    material_with_crystal_sites.analyze()
    return material_with_crystal_sites.get_unique_coordination_numbers(cutoff=cutoff)


def get_coordination_numbers_distribution(
    material: MaterialWithBuildMetadata,
    cutoff: float = 3.0,
) -> Dict[int, int]:
    """
    Get the unique coordination numbers for the provided passivation configuration and cutoff radius.

    Args:
        material (MaterialWithBuildMetadata): The material object.
        cutoff (float): The cutoff radius for defining neighbors.
    Returns:
        set: The unique coordination numbers.
    """
    material_with_crystal_sites = MaterialWithCrystalSites.from_material(material)
    material_with_crystal_sites.analyze()
    coordinatation_numbers = material_with_crystal_sites.get_coordination_numbers(cutoff=cutoff)
    sorted_coordinatation_numbers = sorted(coordinatation_numbers.values)
    distribution = dict(Counter(sorted_coordinatation_numbers))
    return distribution
