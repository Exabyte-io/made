from mat3ra.made.material import Material
from mat3ra.made.tools.analyze import BaseMaterialAnalyzer
from mat3ra.made.tools.convert import from_pymatgen, to_pymatgen

from ..third_party import PymatgenSpacegroupAnalyzer


class LatticeMaterialAnalyzer(BaseMaterialAnalyzer):
    def __init__(self, material: Material):
        super().__init__(material)
        self.spacegroup_analyzer = PymatgenSpacegroupAnalyzer(to_pymatgen(self.material))

    @property
    def get_with_primitive_lattice(self: Material) -> Material:
        """
        Convert a structure to its primitive cell.
        """
        return Material.create(from_pymatgen(self.spacegroup_analyzer.get_primitive_standard_structure()))

    @property
    def get_with_conventional_lattice(self: Material) -> Material:
        """
        Convert a structure to its conventional cell.
        """
        return Material.create(from_pymatgen(self.spacegroup_analyzer.get_conventional_standard_structure()))
