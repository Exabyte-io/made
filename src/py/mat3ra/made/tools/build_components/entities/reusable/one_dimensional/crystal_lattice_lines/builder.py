from typing import Any, Optional, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseSingleBuilder, TypeConfiguration

from ......analyze.lattice_lines import CrystalLatticeLinesMaterialAnalyzer
from ......build.pristine_structures.two_dimensional.slab.utils import get_orthogonal_c_slab
from ......modify import translate_to_z_level
from ......operations.core.unary import supercell
from ..... import MaterialWithBuildMetadata
from .configuration import CrystalLatticeLinesConfiguration


class CrystalLatticeLinesBuilder(BaseSingleBuilder):
    """
    Builder for creating a single crystal lattice line with termination.
    This is similar to CrystalLatticePlanesBuilder but for 1D lines.
    """

    _PostProcessParametersType: Any = None
    use_enforce_convention: bool = True

    def _generate(self, configuration: CrystalLatticeLinesConfiguration) -> Material:
        crystal_lattice_lines_analyzer = CrystalLatticeLinesMaterialAnalyzer(
            material=configuration.crystal, miller_indices_2d=configuration.miller_indices_2d
        )
        miller_supercell_matrix = crystal_lattice_lines_analyzer.miller_supercell_matrix
        miller_supercell_material = supercell(configuration.crystal, miller_supercell_matrix)
        # Lattice returned with vector B being treated as the direction to the surface, needs to be swaped with C
        rotated_material = supercell(miller_supercell_material, [[1, 0, 0], [0, 0, 1], [0, 1, 0]])
        orthogonal_material = get_orthogonal_c_slab(rotated_material)
        return orthogonal_material

    def _enforce_convention(self, material: Union[Material, MaterialWithBuildMetadata]) -> Material:
        if not self.use_enforce_convention:
            return material
        return translate_to_z_level(material, "bottom")

    def _post_process(
        self,
        item: Material,
        post_process_parameters: Optional[Any],
        configuration: Optional[TypeConfiguration] = None,
    ) -> Material:
        item = super()._post_process(item, post_process_parameters, configuration)
        return self._enforce_convention(item)
