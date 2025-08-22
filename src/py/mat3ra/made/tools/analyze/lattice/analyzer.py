from .. import BaseMaterialAnalyzer
from ...build_components.metadata import MaterialWithBuildMetadata
from ...convert import from_pymatgen, to_pymatgen
from ...third_party import PymatgenSpacegroupAnalyzer
from ....lattice import LatticeTypeEnum


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
        try:
            analyzer = PymatgenSpacegroupAnalyzer(
                to_pymatgen(self.material),
                symprec=tolerance,
                angle_tolerance=angle_tolerance,
            )
            lattice_type = analyzer.get_lattice_type()
            spg_symbol = analyzer.get_space_group_symbol()
            
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
