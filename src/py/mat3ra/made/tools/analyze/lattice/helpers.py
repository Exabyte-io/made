from ....lattice import LatticeTypeEnum
from ...build_components.metadata import MaterialWithBuildMetadata
from .analyzer import LatticeMaterialAnalyzer


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


def get_lattice_type(material: MaterialWithBuildMetadata, tolerance=0.2, angle_tolerance=5) -> LatticeTypeEnum:
    """
    Detects the lattice type of the material.

    Args:
        material (MaterialWithBuildMetadata): The material to analyze.
        tolerance (float): Tolerance for lattice parameter comparisons.
        angle_tolerance (float): Tolerance for angle comparisons.

    Returns:
        LatticeTypeEnum: The detected lattice type.
    """
    analyzer = LatticeMaterialAnalyzer(material=material)
    return analyzer.detect_lattice_type(precision=tolerance, angle_tolerance=angle_tolerance)
