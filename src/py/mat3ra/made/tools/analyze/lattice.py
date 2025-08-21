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
    def material_with_conventional_lattice(self: MaterialWithBuildMetadata) -> MaterialWithBuildMetadata:
        """
        Convert a structure to its conventional cell.
        """
        return MaterialWithBuildMetadata.create(
            from_pymatgen(self.spacegroup_analyzer.get_conventional_standard_structure())
        )

    def get_material_with_corrected_orientation(
        self, original_material: MaterialWithBuildMetadata, layer_thickness: float = 1.0
    ) -> MaterialWithBuildMetadata:
        """
        Get the primitive material with corrected orientation if it was flipped.

        Args:
            original_material: The original material before primitivization
            layer_thickness: Thickness of layers for orientation detection

        Returns:
            MaterialWithBuildMetadata: Material with corrected orientation
        """

        basis_analyzer = BasisMaterialAnalyzer(material=self.material)
        if basis_analyzer.is_orientation_flipped(original_material, layer_thickness):
            corrected_material = rotate(self.material, axis=[1, 0, 0], angle=180, rotate_cell=False)
            return corrected_material

        return self.material


def get_material_with_conventional_lattice(material: MaterialWithBuildMetadata) -> MaterialWithBuildMetadata:
    analyzer = LatticeMaterialAnalyzer(material=material)
    return analyzer.material_with_conventional_lattice


def get_material_with_primitive_lattice(
    material: MaterialWithBuildMetadata, return_original_if_not_reduced=False, keep_orientation=True
) -> MaterialWithBuildMetadata:
    analyzer = LatticeMaterialAnalyzer(material=material)
    material_with_primitive_lattice = analyzer.material_with_primitive_lattice
    original_number_of_atoms = material.basis.number_of_atoms
    primitive_structure_number_of_atoms = material_with_primitive_lattice.basis.number_of_atoms

    if original_number_of_atoms == primitive_structure_number_of_atoms:
        if return_original_if_not_reduced:
            return material

    if keep_orientation:
        primitive_analyzer = LatticeMaterialAnalyzer(material=material_with_primitive_lattice)
        material_with_primitive_lattice = primitive_analyzer.get_material_with_corrected_orientation(material)
        print("Orientation corrected after primitive conversion.")

    return material_with_primitive_lattice
