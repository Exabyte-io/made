from collections import Counter
from typing import Union, List, Dict

from mat3ra.made.material import Material
from .configuration import PassivationConfiguration
from .builders import (
    SurfacePassivationBuilder,
    CoordinationBasedPassivationBuilder,
    SurfacePassivationBuilderParameters,
    CoordinationBasedPassivationBuilderParameters,
)
from ...analyze.material import MaterialWithCrystalSites


def create_passivation(
    configuration: PassivationConfiguration,
    builder: Union[SurfacePassivationBuilder, CoordinationBasedPassivationBuilder, None] = None,
) -> Material:
    if builder is None:
        builder = SurfacePassivationBuilder(build_parameters=SurfacePassivationBuilderParameters())
    return builder.get_material(configuration)


def get_unique_coordination_numbers(
    configuration: PassivationConfiguration,
    cutoff: float = 3.0,
) -> List[int]:
    """
    Get the unique coordination numbers for the provided passivation configuration and cutoff radius.

    Args:
        configuration (PassivationConfiguration): The configuration object.
        cutoff (float): The cutoff radius for defining neighbors.
    Returns:
        set: The unique coordination numbers.
    """
    material_with_crystal_sites = MaterialWithCrystalSites.from_material(configuration.slab)
    material_with_crystal_sites.analyze()
    return material_with_crystal_sites.get_unique_coordination_numbers(cutoff=cutoff)


def get_coordination_numbers_distribution(
    configuration: PassivationConfiguration,
    cutoff: float = 3.0,
) -> Dict[int, int]:
    """
    Get the unique coordination numbers for the provided passivation configuration and cutoff radius.

    Args:
        configuration (PassivationConfiguration): The configuration object.
        cutoff (float): The cutoff radius for defining neighbors.
    Returns:
        set: The unique coordination numbers.
    """
    material_with_crystal_sites = MaterialWithCrystalSites.from_material(configuration.slab)
    material_with_crystal_sites.analyze()
    coordinatation_numbers = material_with_crystal_sites.get_coordination_numbers(cutoff=cutoff)
    sorted_coordinatation_numbers = sorted(coordinatation_numbers.values)
    distribution = dict(Counter(sorted_coordinatation_numbers))
    return distribution
