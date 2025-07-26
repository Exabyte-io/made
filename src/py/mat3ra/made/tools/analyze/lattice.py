from mat3ra.made.material import Material
from mat3ra.made.tools.analyze import BaseMaterialAnalyzer
from mat3ra.made.tools.convert import from_pymatgen, to_pymatgen

from ..build import MaterialWithBuildMetadata
from ..third_party import PymatgenSpacegroupAnalyzer


class LatticeMaterialAnalyzer(BaseMaterialAnalyzer):
    @property
    def spacegroup_analyzer(self):
        return PymatgenSpacegroupAnalyzer(to_pymatgen(self.material))

    @property
    def material_with_primitive_lattice(self: Material) -> MaterialWithBuildMetadata:
        """
        Convert a structure to its primitive cell.
        """
        return MaterialWithBuildMetadata.create(
            from_pymatgen(self.spacegroup_analyzer.get_primitive_standard_structure())
        )

    @property
    def material_with_conventional_lattice(self: Material) -> MaterialWithBuildMetadata:
        """
        Convert a structure to its conventional cell.
        """
        return MaterialWithBuildMetadata.create(
            from_pymatgen(self.spacegroup_analyzer.get_conventional_standard_structure())
        )


def get_conventional_material(
    material: MaterialWithBuildMetadata, use_conventional_cell: bool = True
) -> MaterialWithBuildMetadata:
    if use_conventional_cell:
        analyzer = LatticeMaterialAnalyzer(material=material)
        return analyzer.material_with_conventional_lattice
    return material


def get_primitive_material(
    material: MaterialWithBuildMetadata, use_primitive_cell: bool = True
) -> MaterialWithBuildMetadata:
    if use_primitive_cell:
        analyzer = LatticeMaterialAnalyzer(material=material)
        return analyzer.material_with_primitive_lattice
    return material
