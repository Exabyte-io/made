from ....lattice import LatticeTypeEnum
from ...build_components.metadata import MaterialWithBuildMetadata
from ...convert import from_pymatgen, to_pymatgen
from ...third_party import PymatgenSpacegroupAnalyzer
from .. import BaseMaterialAnalyzer
from .utils import detect_lattice_type_from_vectors

PYMATGEN_LATTICE_TYPE_MAP = {
    "cubic": LatticeTypeEnum.CUB,
    "hexagonal": LatticeTypeEnum.HEX,
    "tetragonal": LatticeTypeEnum.TET,
    "rhombohedral": LatticeTypeEnum.RHL,
    "orthorhombic": LatticeTypeEnum.ORC,
    "monoclinic": LatticeTypeEnum.MCL,
    "triclinic": LatticeTypeEnum.TRI,
}


class LatticeMaterialAnalyzer(BaseMaterialAnalyzer):
    @property
    def spacegroup_analyzer(self):
        return PymatgenSpacegroupAnalyzer(to_pymatgen(self.material))

    def detect_lattice_type(self, tolerance=0.1, angle_tolerance=5) -> LatticeTypeEnum:
        """
        Detects the lattice type of the material.

        Args:
            tolerance (float): Tolerance for lattice parameter comparison, in Angstroms.
            angle_tolerance (float): Tolerance for angle comparisons, in degrees.

        Returns:
            LatticeTypeEnum: The detected lattice type.
        """
        # First try vector-based detection for more accurate primitive cell identification
        try:
            vector_based_type = detect_lattice_type_from_vectors(
                self.material.lattice.vector_arrays, 
                tolerance=tolerance, 
                angle_tolerance=angle_tolerance
            )
            
            # If vector-based detection gives a specific result (not TRI), use it
            if vector_based_type != LatticeTypeEnum.TRI:
                return vector_based_type
        except Exception:
            # Fall back to pymatgen if vector-based detection fails
            pass
        
        # Fallback to pymatgen spacegroup analyzer
        lattice_type_str = PymatgenSpacegroupAnalyzer(
            to_pymatgen(self.material),
            symprec=tolerance,
            angle_tolerance=angle_tolerance,
        ).get_lattice_type()

        return PYMATGEN_LATTICE_TYPE_MAP.get(lattice_type_str, LatticeTypeEnum.TRI)

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
