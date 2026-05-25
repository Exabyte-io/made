from typing import Union

from mat3ra.made.material import Material

from ....build.pristine_structures.two_dimensional.slab import SlabConfiguration
from ....build_components.entities.reusable.three_dimensional.supercell.helpers import create_supercell
from ...rdf import RadialDistributionFunction
from .holders import MatchedSubstrateFilmConfigurationHolder


def calculate_interfacial_distance_from_rdf(
    substrate_material: Union[Material, dict, "SlabConfiguration"],
    film_material: Union[Material, dict, "SlabConfiguration"],
    rdf_cutoff: float = 10.0,
    rdf_bin_size: float = 0.1,
    supercell_size: tuple = (3, 3, 3),
) -> float:
    """
    Calculate interfacial distance based on RDF analysis of bulk materials.

    Creates temporary supercells of substrate and film bulk materials,
    calculates their RDFs to find the first peak (nearest neighbor distance),
    and returns the average of these distances as the initial guess for interfacial distance.

    Args:
        substrate_material: Material, dict, or SlabConfiguration for the substrate
        film_material: Material, dict, or SlabConfiguration for the film
        rdf_cutoff: Maximum distance for RDF calculation in Angstroms
        rdf_bin_size: Bin size for RDF histogram in Angstroms
        supercell_size: Size of supercell for RDF analysis (default: 3x3x3)

    Returns:
        float: Calculated interfacial distance in Angstroms
    """

    if isinstance(substrate_material, SlabConfiguration):
        substrate_bulk = substrate_material.atomic_layers.crystal
    elif isinstance(substrate_material, dict):
        substrate_bulk = Material.create(substrate_material)
    else:
        substrate_bulk = substrate_material

    if isinstance(film_material, SlabConfiguration):
        film_bulk = film_material.atomic_layers.crystal
    elif isinstance(film_material, dict):
        film_bulk = Material.create(film_material)
    else:
        film_bulk = film_material

    substrate_supercell = create_supercell(material=substrate_bulk, scaling_factor=list(supercell_size))

    film_supercell = create_supercell(material=film_bulk, scaling_factor=list(supercell_size))

    substrate_rdf = RadialDistributionFunction.from_material(
        substrate_supercell,
        cutoff=rdf_cutoff,
        bin_size=rdf_bin_size,
    )

    film_rdf = RadialDistributionFunction.from_material(
        film_supercell,
        cutoff=rdf_cutoff,
        bin_size=rdf_bin_size,
    )

    substrate_first_peak = substrate_rdf.first_peak_distance
    film_first_peak = film_rdf.first_peak_distance

    interfacial_distance = (substrate_first_peak + film_first_peak) / 2.0

    return interfacial_distance


__all__ = [
    "MatchedSubstrateFilmConfigurationHolder",
    "calculate_interfacial_distance_from_rdf",
]
