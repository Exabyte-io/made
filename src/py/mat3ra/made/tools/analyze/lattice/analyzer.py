from ....lattice import LatticeTypeEnum
from ...build_components.metadata import MaterialWithBuildMetadata
from ...convert import from_pymatgen, to_pymatgen
from ...operations.core.unary import rotate
from ...third_party import PymatgenSpacegroupAnalyzer
from .. import BaseMaterialAnalyzer
from ..basis import BasisMaterialAnalyzer


class LatticeMaterialAnalyzer(BaseMaterialAnalyzer):
    precision: float = 0.1
    angle_tolerance: float = 5.0

    @property
    def spacegroup_analyzer(self):
        return PymatgenSpacegroupAnalyzer(
            to_pymatgen(self.material), symprec=self.precision, angle_tolerance=self.angle_tolerance
        )

    def detect_lattice_type(self, precision=0.1, angle_tolerance=5) -> LatticeTypeEnum:
        """
        Detects the lattice type of the material.

        Args:
            precision (float): Tolerance for lattice parameter comparison, in Angstroms.
            angle_tolerance (float): Tolerance for angle comparisons, in degrees.

        Returns:
            LatticeTypeEnum: The detected lattice type.
        """
        self.precision = precision
        self.angle_tolerance = angle_tolerance
        try:
            lattice_type = self.spacegroup_analyzer.get_lattice_type()
            spg_symbol = self.spacegroup_analyzer.get_space_group_symbol()

            # Enhanced detection using space group symbol
            if lattice_type == "cubic":
                if "P" in spg_symbol:
                    return LatticeTypeEnum.CUB
                elif "F" in spg_symbol:
                    return LatticeTypeEnum.FCC
                elif "I" in spg_symbol:
                    return LatticeTypeEnum.BCC
            elif lattice_type == "tetragonal":
                if "P" in spg_symbol:
                    return LatticeTypeEnum.TET
                elif "I" in spg_symbol:
                    return LatticeTypeEnum.BCT
            elif lattice_type == "orthorhombic":
                if "P" in spg_symbol:
                    return LatticeTypeEnum.ORC
                elif "F" in spg_symbol:
                    return LatticeTypeEnum.ORCF
                elif "I" in spg_symbol:
                    return LatticeTypeEnum.ORCI
                elif "C" in spg_symbol:
                    return LatticeTypeEnum.ORCC
            elif lattice_type == "hexagonal":
                return LatticeTypeEnum.HEX
            elif lattice_type == "rhombohedral":
                return LatticeTypeEnum.RHL
            elif lattice_type == "monoclinic":
                if "P" in spg_symbol:
                    return LatticeTypeEnum.MCL
                elif "C" in spg_symbol:
                    return LatticeTypeEnum.MCLC

        except Exception:
            return LatticeTypeEnum.TRI

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
