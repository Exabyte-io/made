from typing import List

from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer

from ...material import Material
from ...utils import get_center_of_coordinates
from ..analyze import BaseMaterialAnalyzer
from ..analyze.coordination import get_voronoi_nearest_neighbors_atom_indices
from ..analyze.other import get_closest_site_id_from_coordinate
from ..build.defect.slab.helpers import recreate_slab_with_fractional_layers
from ..build.supercell import create_supercell
from ..convert import to_pymatgen
from ..modify import filter_by_condition_on_coordinates
from ..third_party import PymatgenVoronoiInterstitialGenerator
from ..utils import get_distance_between_coordinates, transform_coordinate_to_supercell


class CrystalSiteAnalyzer(BaseMaterialAnalyzer):
    coordinate: List[float]

    @property
    def exact_coordinate(self) -> List[float]:
        return self.coordinate

    @property
    def closest_site_coordinate(self) -> List[float]:
        site_id = get_closest_site_id_from_coordinate(self.material, self.coordinate)
        return self.material.coordinates_array[site_id]

    @property
    def equidistant_coordinate(self) -> List[float]:
        """
        Compute a coordinate that is equidistant from the nearest atoms to the target coordinate. Useful for adatom.

        This method works by:
        1. Creating a 3x3x1 supercell of the original material to include atoms in PBC.
        2. Transforming the target coordinate into the supercell's coordinate system.
        3. Using Voronoi tessellation to find the indices of the nearest atoms in the supercell.
        4. Taking the geometric center (average) of these neighboring atoms' coordinates as the equidistant point.
        5. Setting the z-coordinate to the original value (useful for 2D/slab systems).
        6. Transforming the equidistant coordinate back to the original cell's coordinate system.

        The [3, 3, 1] supercell ensures robust neighbor search in x and y, but not in z (for 2D/slab systems).
        """
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

        # Filter out atoms that are too close to the z boundaries of the supercell
        supercell_material = filter_by_condition_on_coordinates(supercell_material, lambda c: 1e-2 < c[2] < 1 - 1e-2)
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
    """
    Attributes:
       clustering_tol: Tolerance for clustering the Voronoi nodes.
       min_dist: Minimum distance between an interstitial and the nearest atom.
       ltol: Tolerance for lattice matching.
       stol: Tolerance for structure matching.
       angle_tol: Angle tolerance for structure matching.
    """

    clustering_tol: float = 0.5
    min_dist: float = 0.9
    ltol: float = 0.2
    stol: float = 0.3
    angle_tol: float = 5

    @property
    def voronoi_site_coordinate(self) -> List[float]:
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


class AdatomCrystalSiteAnalyzer(CrystalSiteAnalyzer, SlabMaterialAnalyzer):
    """
    This analyzer is used to find a new crystal site closest to the given coordinate.
    By creating the slab with additional layers and finding the closest site in the new slab.
    """

    @property
    def slab_analyzer(self) -> SlabMaterialAnalyzer:
        return SlabMaterialAnalyzer(material=self.material)

    @property
    def new_crystal_site_coordinate(self) -> List[float]:
        new_slab = self.get_additional_layers_slab()
        analyzer = CrystalSiteAnalyzer(material=new_slab, coordinate=self.coordinate)
        return analyzer.closest_site_coordinate

    def get_additional_layers_slab(self) -> Material:
        # number of layers can be calculated to allow for coordinate to resolved there
        slab_with_one_layer = recreate_slab_with_fractional_layers(self.slab_analyzer.material, number_of_layers=1)
        return slab_with_one_layer
