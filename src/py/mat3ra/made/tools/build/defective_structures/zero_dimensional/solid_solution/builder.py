import random
from typing import List, Type, Union

import numpy as np
from mat3ra.made.material import Material

from .....analyze.other import get_chemical_formula_empirical
from .....build_components import (
    BaseBuilderParameters,
    BaseSingleBuilder,
    MaterialWithBuildMetadata,
    TypeConfiguration,
)
from .....build_components.entities.reusable.three_dimensional.supercell.helpers import create_supercell
from .configuration import SolidSolutionConfiguration


def _get_cartesian_coordinates(material: Union[Material, MaterialWithBuildMetadata]) -> np.ndarray:
    cloned = material.clone()
    cloned.to_cartesian()
    return np.array(cloned.basis.coordinates.values)


def _minimum_image_distances(
    frac_coords: np.ndarray,
    lattice_vectors: np.ndarray,
) -> np.ndarray:
    """
    Compute pairwise minimum image distances for fractional coordinates under PBC.

    Args:
        frac_coords (np.ndarray): Fractional coordinates, shape (N, 3).
        lattice_vectors (np.ndarray): Lattice vectors as rows, shape (3, 3).

    Returns:
        np.ndarray: Symmetric distance matrix, shape (N, N).
    """
    n = len(frac_coords)
    distances = np.zeros((n, n))
    for i in range(n):
        delta = frac_coords[i] - frac_coords[i + 1 :]
        delta -= np.round(delta)
        cart = delta @ lattice_vectors
        dists = np.linalg.norm(cart, axis=1)
        distances[i, i + 1 :] = dists
        distances[i + 1 :, i] = dists
    return distances


def _farthest_point_sampling(
    source_indices: List[int],
    frac_coords: np.ndarray,
    lattice_vectors: np.ndarray,
    n_select: int,
    seed: int = None,
) -> List[int]:
    """
    Select n_select sites from source_indices using Farthest Point Sampling.

    Iteratively picks the site whose minimum distance to all already-selected
    sites is largest, producing a maximally dispersed subset under PBC.
    """
    if n_select <= 0:
        return []
    if n_select >= len(source_indices):
        return sorted(source_indices)

    sublattice_coords = frac_coords[source_indices]
    dist_matrix = _minimum_image_distances(sublattice_coords, lattice_vectors)
    n = len(source_indices)

    rng = random.Random(seed)
    start = rng.randrange(n)
    selected = [start]
    min_dist = dist_matrix[start].copy()
    min_dist[start] = -1.0

    for _ in range(n_select - 1):
        next_site = int(np.argmax(min_dist))
        selected.append(next_site)
        np.minimum(min_dist, dist_matrix[next_site], out=min_dist)
        min_dist[next_site] = -1.0

    return sorted([source_indices[s] for s in selected])


class SolidSolutionBuilderParameters(BaseBuilderParameters):
    site_selection_method: str = "random"


class SolidSolutionBuilder(BaseSingleBuilder):
    _ConfigurationType: Type[SolidSolutionConfiguration] = SolidSolutionConfiguration
    _BuildParametersType: Type[SolidSolutionBuilderParameters] = SolidSolutionBuilderParameters
    _DefaultBuildParameters: SolidSolutionBuilderParameters = SolidSolutionBuilderParameters()

    def _generate(self, configuration: SolidSolutionConfiguration) -> MaterialWithBuildMetadata:
        supercell = create_supercell(
            material=configuration.crystal,
            scaling_factor=configuration.supercell_dimensions,
        )
        material = MaterialWithBuildMetadata.create(supercell.to_dict())

        elements = material.basis.elements.values
        source_indices = [i for i, el in enumerate(elements) if el == configuration.source_element]

        n_replace = round(configuration.concentration * len(source_indices))
        n_replace = max(0, min(n_replace, len(source_indices)))

        method = configuration.site_selection_method

        if method == "uniform":
            material.to_crystal()
            frac_coords = np.array(material.basis.coordinates.values)
            lattice_vectors = np.array(material.lattice.vector_arrays)
            selected_indices = _farthest_point_sampling(
                source_indices, frac_coords, lattice_vectors, n_replace, seed=configuration.seed
            )
        else:
            rng = random.Random(configuration.seed)
            selected_indices = sorted(rng.sample(source_indices, n_replace))

        for idx in selected_indices:
            material.basis.elements.values[idx] = configuration.target_element

        return material

    def _update_material_name(
        self, material: Union[Material, MaterialWithBuildMetadata], configuration: TypeConfiguration
    ) -> MaterialWithBuildMetadata:
        formula = get_chemical_formula_empirical(material)
        material.name = f"{formula} - Solid Solution"
        return material
