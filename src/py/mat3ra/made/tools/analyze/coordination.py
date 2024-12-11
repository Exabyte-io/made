from typing import List, Optional

from mat3ra.made.material import Material

from ..convert import to_pymatgen
from ..third_party import PymatgenVoronoiNN


def get_voronoi_nearest_neighbors_atom_indices(
    material: Material,
    coordinate: Optional[List[float]] = None,
    tolerance: float = 0.1,
    cutoff: float = 13.0,
) -> Optional[List[int]]:
    """
    Returns the indices of direct neighboring atoms to a specified position in the material using Voronoi tessellation.

    Args:
        material (Material): The material object to find neighbors in.
        coordinate (List[float]): The position to find neighbors for.
        tolerance (float): tolerance parameter for near-neighbor finding. Faces that are smaller than tol fraction
            of the largest face are not included in the tessellation. (default: 0.1).
            as per: https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.local_env.VoronoiNN
        cutoff (float): The cutoff radius for identifying neighbors, in angstroms.

    Returns:
        List[int]: A list of indices of neighboring atoms, or an empty list if no neighbors are found.
    """
    if coordinate is None:
        coordinate = [0, 0, 0]
    structure = to_pymatgen(material)
    voronoi_nn = PymatgenVoronoiNN(
        tol=tolerance,
        cutoff=cutoff,
        weight="solid_angle",
        extra_nn_info=False,
        compute_adj_neighbors=True,
    )
    coordinates = material.basis.coordinates
    site_index = coordinates.get_element_id_by_value(coordinate)
    remove_dummy_atom = False
    if site_index is None:
        structure.append("X", coordinate, validate_proximity=False)
        site_index = len(structure.sites) - 1
        remove_dummy_atom = True
    try:
        neighbors = voronoi_nn.get_nn_info(structure, site_index)
    except ValueError:
        return None
    neighboring_atoms_pymatgen_ids = [n["site_index"] for n in neighbors]
    if remove_dummy_atom:
        structure.remove_sites([-1])

    all_coordinates = material.basis.coordinates
    all_coordinates.filter_by_indices(neighboring_atoms_pymatgen_ids)
    return all_coordinates.ids
