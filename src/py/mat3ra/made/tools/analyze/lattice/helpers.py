from ....lattice import LatticeTypeEnum
from ...build_components.metadata import MaterialWithBuildMetadata
from .analyzer import LatticeMaterialAnalyzer


def get_material_with_conventional_lattice(material: MaterialWithBuildMetadata) -> MaterialWithBuildMetadata:
    analyzer = LatticeMaterialAnalyzer(material=material)
    return analyzer.material_with_conventional_lattice


def get_material_with_primitive_lattice(
    material: MaterialWithBuildMetadata, return_original_if_not_reduced=False, keep_orientation=True
) -> MaterialWithBuildMetadata:
    analyzer = LatticeMaterialAnalyzer(material=material)
    return analyzer.get_material_with_primitive_lattice_standard(return_original_if_not_reduced, keep_orientation)


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
