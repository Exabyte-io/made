from typing import List

from mat3ra.made.tools.analyze import BaseMaterialAnalyzer
from mat3ra.made.tools.analyze.coordination import get_voronoi_nearest_neighbors_atom_indices
from mat3ra.made.tools.analyze.other import get_closest_site_id_from_coordinate
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.convert import to_pymatgen
from mat3ra.made.tools.third_party import PymatgenVoronoiInterstitialGenerator
from mat3ra.made.tools.utils import get_distance_between_coordinates, transform_coordinate_to_supercell
from mat3ra.made.utils import get_center_of_coordinates


class CrystalSiteAnalyzer(BaseMaterialAnalyzer):
    coordinate: List[float]

    @property
    def coordinate_resolution(self) -> List[float]:
        return self.coordinate

    @property
    def closest_site_resolution(self) -> List[float]:
        site_id = get_closest_site_id_from_coordinate(self.material, self.coordinate)
        return self.material.coordinates_array[site_id]

    @property
    def new_crystal_site_resolution(self) -> List[float]:
        return self.coordinate

    @property
    def equidistant_resolution(self) -> List[float]:
        scaling_factor = [3, 3, 1]
        translation_vector = [1 / 3, 1 / 3, 0]
        supercell_material = create_supercell(self.material, scaling_factor=scaling_factor)

        coordinate_in_supercell = transform_coordinate_to_supercell(
            coordinate=self.coordinate, scaling_factor=scaling_factor, translation_vector=translation_vector
        )

        neighboring_atoms_ids_in_supercell = get_voronoi_nearest_neighbors_atom_indices(
            material=supercell_material, coordinate=coordinate_in_supercell
        )

        if neighboring_atoms_ids_in_supercell is None:
            raise ValueError("No neighboring atoms found for equidistant calculation.")

        isolated_neighboring_atoms_basis = supercell_material.basis.model_copy()
        isolated_neighboring_atoms_basis.coordinates.filter_by_ids(neighboring_atoms_ids_in_supercell)
        equidistant_coordinate_in_supercell = get_center_of_coordinates(
            isolated_neighboring_atoms_basis.coordinates.values
        )
        equidistant_coordinate_in_supercell[2] = self.coordinate[2]

        return transform_coordinate_to_supercell(
            equidistant_coordinate_in_supercell, scaling_factor, translation_vector, reverse=True
        )


class VoronoiCrystalSiteAnalyzer(CrystalSiteAnalyzer):
    clustering_tol: float = 0.5
    min_dist: float = 0.9
    ltol: float = 0.2
    stol: float = 0.3
    angle_tol: float = 5

    @property
    def voronoi_site_resolution(self) -> List[float]:
        """Voronoi site resolution using configured parameters."""
        return self._get_voronoi_site_resolution()

    def _get_voronoi_site_resolution(self) -> List[float]:
        pymatgen_structure = to_pymatgen(self.material)

        voronoi_gen = PymatgenVoronoiInterstitialGenerator(
            clustering_tol=self.clustering_tol,
            min_dist=self.min_dist,
            ltol=self.ltol,
            stol=self.stol,
            angle_tol=self.angle_tol,
        )

        interstitials = list(voronoi_gen.generate(structure=pymatgen_structure, insert_species=["Si"]))

        if not interstitials:
            raise ValueError("No Voronoi interstitial sites found.")

        closest_interstitial = min(
            interstitials,
            key=lambda interstitial: get_distance_between_coordinates(interstitial.site.frac_coords, self.coordinate),
        )

        return closest_interstitial.site.frac_coords.tolist()
