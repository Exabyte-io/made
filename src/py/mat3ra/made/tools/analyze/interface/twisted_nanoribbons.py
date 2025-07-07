from typing import List, Tuple

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface import InterfaceAnalyzer
from mat3ra.made.tools.analyze.interface.utils.holders import MatchedSubstrateFilmConfigurationHolder
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from mat3ra.made.tools.modify import rotate, translate_to_center


class TwistedNanoribbonsInterfaceAnalyzer(InterfaceAnalyzer):
    """
    Analyzer for creating twisted interfaces between two nanoribbons.

    Takes two nanoribbons (or any 0D materials) and returns them with:
    - Same lattice vectors
    - One rotated by the specified angle
    - Nanoribbons centered
    - Gap applied using SlabStrainedSupercellWithGapConfiguration

    Args:

        angle (float): Twist angle in degrees.
        vacuum_x (float): Vacuum along x on both sides, in Angstroms.
        vacuum_y (float): Vacuum along y on both sides, in Angstroms.
    """

    angle: float = 0.0
    vacuum_x: float = 5.0
    vacuum_y: float = 5.0

    @property
    def nanoribbon1(self):
        return self.substrate_slab_configuration.atomic_layers.crystal

    @property
    def nanoribbon2(self):
        return self.film_slab_configuration.atomic_layers.crystal

    def _create_primitive_slab_from_nanoribbon(self, nanoribbon: Material) -> SlabConfiguration:
        return SlabConfiguration.from_parameters(
            material_or_dict=nanoribbon,
            miller_indices=(0, 0, 1),
            number_of_layers=1,
            vacuum=0.0,
        )

    def _match_lattice_vectors(self, material1: Material, material2: Material) -> Tuple[Material, Material]:
        lattice1 = material1.lattice.vector_arrays
        lattice2 = material2.lattice.vector_arrays

        n_atoms_1 = len(material1.basis.elements.values)
        n_atoms_2 = len(material2.basis.elements.values)
        if n_atoms_1 >= n_atoms_2:
            target_lattice = lattice1
            material1_matched = material1
            material2_matched = material2.clone()
            material2_matched.set_lattice(material1.lattice)
        else:
            target_lattice = lattice2
            material1_matched = material1.clone()
            material1_matched.set_lattice(material2.lattice)
            material2_matched = material2

        return material1_matched, material2_matched

    def get_strained_configuration(self) -> MatchedSubstrateFilmConfigurationHolder:
        """
        Get strained configurations for the twisted nanoribbon interface.

        Returns:
            List[MatchedSubstrateFilmConfigurationHolder]: List of strained configurations
        """
        centered_nanoribbon1 = translate_to_center(self.nanoribbon1, axes=["x", "y"])
        centered_nanoribbon2 = translate_to_center(self.nanoribbon2, axes=["x", "y"])

        matched_nanoribbon1, matched_nanoribbon2 = self._match_lattice_vectors(
            centered_nanoribbon1, centered_nanoribbon2
        )

        rotated_nanoribbon2 = rotate(matched_nanoribbon2, [0, 0, 1], self.angle, wrap=False)

        substrate_slab_config = self._create_primitive_slab_from_nanoribbon(matched_nanoribbon1)
        film_slab_config = self._create_primitive_slab_from_nanoribbon(rotated_nanoribbon2)

        return self.create_matched_configuration_holder(substrate_slab_config, film_slab_config, match_id=0)

    @property
    def substrate_nanoribbon_configuration(self):
        return self.get_strained_configuration().substrate_configuration

    @property
    def film_nanoribbon_configuration(self):
        return self.get_strained_configuration().film_configuration
