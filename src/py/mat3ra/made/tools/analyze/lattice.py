from ..build_components.metadata import MaterialWithBuildMetadata
from ..convert import from_pymatgen, to_pymatgen
from ..operations.core.unary import rotate
from ..third_party import PymatgenSpacegroupAnalyzer
from . import BaseMaterialAnalyzer
from .basis import BasisMaterialAnalyzer


class LatticeMaterialAnalyzer(BaseMaterialAnalyzer):
    @property
    def spacegroup_analyzer(self):
        return PymatgenSpacegroupAnalyzer(to_pymatgen(self.material))

    @property
    def material_with_primitive_lattice(self: MaterialWithBuildMetadata) -> MaterialWithBuildMetadata:
        """
        Convert a structure to its primitive cell.
        """
        return MaterialWithBuildMetadata.create(
            from_pymatgen(self.spacegroup_analyzer.get_primitive_standard_structure())
        )

    @property
    def material_with_primitive_lattice_standard(self) -> MaterialWithBuildMetadata:
        """
        Get material with primitive lattice and standard orientation correction.
        Uses default parameters: return_original_if_not_reduced=False, keep_orientation=True
        """
        return self.get_material_with_primitive_lattice_standard()

    @property
    def material_with_conventional_lattice(self: MaterialWithBuildMetadata) -> MaterialWithBuildMetadata:
        """
        Convert a structure to its conventional cell.
        """
        return MaterialWithBuildMetadata.create(
            from_pymatgen(self.spacegroup_analyzer.get_conventional_standard_structure())
        )

    def get_material_with_primitive_lattice_standard(
        self, return_original_if_not_reduced: bool = False, keep_orientation: bool = True, layer_thickness: float = 1.0
    ) -> MaterialWithBuildMetadata:
        """
        Get material with primitive lattice and optional orientation correction to be standardized.

        Args:
            return_original_if_not_reduced: If True, return original material when no reduction occurs
            keep_orientation: If True, correct orientation after primitive conversion
            layer_thickness: Thickness of layers for orientation detection

        Returns:
            MaterialWithBuildMetadata: Material with primitive lattice
        """
        material_with_primitive_lattice = self.material_with_primitive_lattice
        original_number_of_atoms = self.material.basis.number_of_atoms
        primitive_structure_number_of_atoms = material_with_primitive_lattice.basis.number_of_atoms

        if original_number_of_atoms == primitive_structure_number_of_atoms:
            if return_original_if_not_reduced:
                return self.material

        if keep_orientation:
            basis_analyzer = BasisMaterialAnalyzer(material=material_with_primitive_lattice)
            if basis_analyzer.is_orientation_flipped(self.material, layer_thickness):
                material_with_primitive_lattice = rotate(
                    material_with_primitive_lattice, axis=[1, 0, 0], angle=180, rotate_cell=False
                )
                print("Orientation corrected after primitive conversion.")

        return material_with_primitive_lattice


def get_material_with_conventional_lattice(material: MaterialWithBuildMetadata) -> MaterialWithBuildMetadata:
    analyzer = LatticeMaterialAnalyzer(material=material)
    return analyzer.material_with_conventional_lattice


def get_material_with_primitive_lattice(
    material: MaterialWithBuildMetadata, return_original_if_not_reduced=False, keep_orientation=True
) -> MaterialWithBuildMetadata:
    analyzer = LatticeMaterialAnalyzer(material=material)
    return analyzer.get_material_with_primitive_lattice_standard(return_original_if_not_reduced, keep_orientation)
