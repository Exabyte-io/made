import math

import numpy as np
from mat3ra.code.array_with_ids import ArrayWithIds

from .. import BaseMaterialAnalyzer
from ..lattice_swap_analyzer import MaterialLatticeSwapAnalyzer
from ...build_components.metadata import MaterialWithBuildMetadata
from ...convert import from_pymatgen, to_pymatgen
from ...third_party import PymatgenSpacegroupAnalyzer
from ....lattice import LatticeTypeEnum


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
        self,
        return_original_if_not_reduced: bool = False,
        keep_orientation: bool = True,
        layer_thickness: float = 1.0,
        rotation_detection_threshold: float = 0.05,
    ) -> MaterialWithBuildMetadata:
        """
        Get material with primitive lattice standardized according to IUCr conventions.

        Args:
            return_original_if_not_reduced: If True, return original material when no reduction occurs
            keep_orientation: If True, detect and reverse lattice parameter swaps to preserve original orientation
            layer_thickness: Unused (kept for compatibility)
            rotation_detection_threshold: Unused (kept for compatibility)

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
            swap_analyzer = MaterialLatticeSwapAnalyzer(material=material_with_primitive_lattice)
            material_with_primitive_lattice = swap_analyzer.get_corrected_material(
                self.material, layer_thickness=layer_thickness, threshold=rotation_detection_threshold
            )

        material_with_primitive_lattice.metadata = self.material.metadata

        return material_with_primitive_lattice

    def _canonicalize_in_plane_primitive_material(
        self, material: MaterialWithBuildMetadata, tolerance: float = 1e-8
    ) -> MaterialWithBuildMetadata:
        if not self._should_canonicalize_in_plane(material, tolerance):
            return material

        canonical_material = material.clone()
        original_is_cartesian = canonical_material.basis.is_in_cartesian_units
        canonical_material.to_crystal()

        coordinates = np.array(canonical_material.coordinates_array, dtype=float)
        if len(coordinates) == 0:
            return canonical_material

        elements = list(canonical_material.basis.elements.values)
        labels = list(canonical_material.basis.labels.values)
        constraints = list(canonical_material.basis.constraints.values)

        best_representation = None
        best_coordinates = None
        best_order = None

        for operation in self._get_in_plane_symmetry_operations(material, tolerance):
            transformed_coordinates = coordinates.copy()
            for index, coordinate in enumerate(coordinates):
                x_value, y_value = operation(coordinate[0], coordinate[1])
                transformed_coordinates[index, 0] = x_value
                transformed_coordinates[index, 1] = y_value

            x_candidates = sorted(
                {
                    self._normalize_fractional_coordinate(coord[0], tolerance)
                    for coord in transformed_coordinates
                }
            )
            y_candidates = sorted(
                {
                    self._normalize_fractional_coordinate(coord[1], tolerance)
                    for coord in transformed_coordinates
                }
            )

            for x_shift in x_candidates:
                for y_shift in y_candidates:
                    shifted_coordinates = transformed_coordinates.copy()
                    shifted_coordinates[:, 0] = [
                        self._normalize_fractional_coordinate(value - x_shift, tolerance)
                        for value in shifted_coordinates[:, 0]
                    ]
                    shifted_coordinates[:, 1] = [
                        self._normalize_fractional_coordinate(value - y_shift, tolerance)
                        for value in shifted_coordinates[:, 1]
                    ]

                    order = sorted(
                        range(len(shifted_coordinates)),
                        key=lambda index: self._canonical_atom_key(
                            elements[index],
                            labels[index] if labels else None,
                            shifted_coordinates[index],
                            tolerance,
                        ),
                    )
                    representation = tuple(
                        self._canonical_atom_key(
                            elements[index],
                            labels[index] if labels else None,
                            shifted_coordinates[index],
                            tolerance,
                        )
                        for index in order
                    )

                    if best_representation is None or representation < best_representation:
                        best_representation = representation
                        best_coordinates = shifted_coordinates
                        best_order = order

        if best_coordinates is None or best_order is None:
            return material

        canonical_material.basis.coordinates.values = best_coordinates[best_order].tolist()
        canonical_material.basis.elements = ArrayWithIds.from_values(
            [elements[index] for index in best_order]
        )
        if labels:
            canonical_material.basis.labels = ArrayWithIds.from_values(
                [labels[index] for index in best_order]
            )
        if constraints and len(constraints) == len(best_order):
            canonical_material.basis.constraints = ArrayWithIds.from_values(
                [constraints[index] for index in best_order]
            )

        if original_is_cartesian:
            canonical_material.to_cartesian()

        return canonical_material

    def _should_canonicalize_in_plane(self, material: MaterialWithBuildMetadata, tolerance: float) -> bool:
        return (
            math.isclose(material.lattice.alpha, 90.0, abs_tol=tolerance)
            and math.isclose(material.lattice.beta, 90.0, abs_tol=tolerance)
            and math.isclose(material.lattice.gamma, 90.0, abs_tol=tolerance)
            and math.isclose(material.lattice.a, material.lattice.b, rel_tol=0.0, abs_tol=tolerance)
            and material.lattice.c > max(material.lattice.a, material.lattice.b)
        )

    def _get_in_plane_symmetry_operations(self, material: MaterialWithBuildMetadata, tolerance: float):
        operations = [
            lambda x_value, y_value: (x_value, y_value),
            lambda x_value, y_value: (-x_value, y_value),
            lambda x_value, y_value: (x_value, -y_value),
            lambda x_value, y_value: (-x_value, -y_value),
        ]

        if math.isclose(material.lattice.a, material.lattice.b, rel_tol=0.0, abs_tol=tolerance):
            operations.extend(
                [
                    lambda x_value, y_value: (y_value, x_value),
                    lambda x_value, y_value: (-y_value, x_value),
                    lambda x_value, y_value: (y_value, -x_value),
                    lambda x_value, y_value: (-y_value, -x_value),
                ]
            )

        return operations

    @staticmethod
    def _normalize_fractional_coordinate(value: float, tolerance: float) -> float:
        normalized_value = float(value % 1.0)
        if abs(normalized_value) < tolerance or abs(normalized_value - 1.0) < tolerance:
            return 0.0
        return normalized_value

    def _canonical_atom_key(self, element, label, coordinate, tolerance: float):
        rounded_coordinate = tuple(
            round(self._normalize_fractional_coordinate(value, tolerance), 8)
            for value in coordinate
        )
        label_key = "" if label is None else str(label)
        return (element, label_key, rounded_coordinate[2], rounded_coordinate[0], rounded_coordinate[1])
