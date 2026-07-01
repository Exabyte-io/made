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


def _minimum_image_distances(frac_coords: np.ndarray, lattice_vectors: np.ndarray) -> np.ndarray:
    n = len(frac_coords)
    distances = np.zeros((n, n))
    for i in range(n):
        delta = frac_coords[i] - frac_coords[i + 1 :]
        delta -= np.round(delta)
        dists = np.linalg.norm(delta @ lattice_vectors, axis=1)
        distances[i, i + 1 :] = dists
        distances[i + 1 :, i] = dists
    return distances


def _select_sites_uniform(
    material: Union[Material, MaterialWithBuildMetadata],
    source_indices: List[int],
    n_select: int,
    seed: int = None,
) -> List[int]:
    """
    Select n_select sites from source_indices using Farthest Point Sampling.

    Iteratively picks the site whose minimum distance to all already-selected
    sites is largest, producing a maximally dispersed subset under PBC.

    Args:
        material: Material with lattice and coordinates.
        source_indices: Indices of candidate sites in the material basis.
        n_select: Number of sites to select.
        seed: Random seed for the initial site choice.

    Returns:
        Sorted list of selected site indices.
    """
    if n_select <= 0:
        return []
    if n_select >= len(source_indices):
        return sorted(source_indices)

    material.to_crystal()
    frac_coords = np.array(material.basis.coordinates.values)
    lattice_vectors = np.array(material.lattice.vector_arrays)
    dist_matrix = _minimum_image_distances(frac_coords[source_indices], lattice_vectors)

    rng = random.Random(seed)
    start = rng.randrange(len(source_indices))
    selected = [start]
    min_dist = dist_matrix[start].copy()
    min_dist[start] = -1.0

    for _ in range(n_select - 1):
        next_idx = int(np.argmax(min_dist))
        selected.append(next_idx)
        np.minimum(min_dist, dist_matrix[next_idx], out=min_dist)
        min_dist[next_idx] = -1.0

    return sorted([source_indices[i] for i in selected])


class SolidSolutionBuilderParameters(BaseBuilderParameters):
    site_selection_method: str = "uniform"


class SolidSolutionBuilder(BaseSingleBuilder):
    _ConfigurationType: Type[SolidSolutionConfiguration] = SolidSolutionConfiguration
    _BuildParametersType: Type[SolidSolutionBuilderParameters] = SolidSolutionBuilderParameters
    _DefaultBuildParameters: SolidSolutionBuilderParameters = SolidSolutionBuilderParameters()

    def _generate(self, configuration: SolidSolutionConfiguration) -> MaterialWithBuildMetadata:
        supercell = create_supercell(material=configuration.crystal, scaling_factor=configuration.supercell_dimensions)
        material = MaterialWithBuildMetadata.create(supercell.to_dict())

        elements = material.basis.elements.values
        source_indices = [i for i, el in enumerate(elements) if el == configuration.source_element]
        n_replace = max(0, min(round(configuration.concentration * len(source_indices)), len(source_indices)))

        if configuration.site_selection_method == "uniform":
            selected = _select_sites_uniform(material, source_indices, n_replace, configuration.seed)
        else:
            rng = random.Random(configuration.seed)
            selected = sorted(rng.sample(source_indices, n_replace))

        for idx in selected:
            elements[idx] = configuration.target_element

        return material

    def _update_material_name(
        self, material: Union[Material, MaterialWithBuildMetadata], configuration: TypeConfiguration
    ) -> MaterialWithBuildMetadata:
        formula = get_chemical_formula_empirical(material)
        material.name = f"{formula} - Solid Solution"
        return material
