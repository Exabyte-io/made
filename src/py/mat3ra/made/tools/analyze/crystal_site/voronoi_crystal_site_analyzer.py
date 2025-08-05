from typing import List

from ...convert import to_pymatgen
from ...third_party import PymatgenVoronoiInterstitialGenerator
from ...utils import get_distance_between_coordinates
from .crystal_site_analyzer import CrystalSiteAnalyzer


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
