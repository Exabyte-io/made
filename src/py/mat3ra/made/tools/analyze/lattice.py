from ..build_components.metadata import MaterialWithBuildMetadata
from ..convert import from_pymatgen, to_pymatgen
from ..third_party import PymatgenSpacegroupAnalyzer
from . import BaseMaterialAnalyzer


class LatticeMaterialAnalyzer(BaseMaterialAnalyzer):
    @property
    def spacegroup_analyzer(self):
        return PymatgenSpacegroupAnalyzer(to_pymatgen(self.material))

    def detect_lattice_type(self, tolerance=0.2, angle_tolerance=5):
        """
        Detects the lattice type of the material.

        Args:
            tolerance (float): Tolerance for lattice parameter comparisons.
            angle_tolerance (float): Tolerance for angle comparisons.

        Returns:
            str: The detected lattice type.
        """
        return PymatgenSpacegroupAnalyzer(
            to_pymatgen(self.material),
            symprec=tolerance,
            angle_tolerance=angle_tolerance,
        ).get_lattice_type()

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


def get_material_with_conventional_lattice(material: MaterialWithBuildMetadata) -> MaterialWithBuildMetadata:
    analyzer = LatticeMaterialAnalyzer(material=material)
    return analyzer.material_with_conventional_lattice


def get_material_with_primitive_lattice(
    material: MaterialWithBuildMetadata, return_original_if_not_reduced=False
) -> MaterialWithBuildMetadata:
    analyzer = LatticeMaterialAnalyzer(material=material)
    material_with_primitive_lattice = analyzer.material_with_primitive_lattice
    original_number_of_atoms = material.basis.number_of_atoms
    primitive_structure_number_of_atoms = material_with_primitive_lattice.basis.number_of_atoms
    if original_number_of_atoms == primitive_structure_number_of_atoms:
        # Not reduced, return original material if requested, to avoid unnecessary editions
        if return_original_if_not_reduced:
            return material
    # Reduced, return the primitive structure
    return material_with_primitive_lattice


def get_lattice_type(material: MaterialWithBuildMetadata, tolerance=0.2, angle_tolerance=5) -> str:
    """
    Detects the lattice type of the material.

    Args:
        material (MaterialWithBuildMetadata): The material to analyze.
        tolerance (float): Tolerance for lattice parameter comparisons.
        angle_tolerance (float): Tolerance for angle comparisons.

    Returns:
        str: The detected lattice type.
    """
    analyzer = LatticeMaterialAnalyzer(material=material)
    return analyzer.detect_lattice_type(tolerance=tolerance, angle_tolerance=angle_tolerance)
