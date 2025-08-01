import numpy as np
from mat3ra.made.tools.analyze import BaseMaterialAnalyzer
from mat3ra.made.tools.convert import from_pymatgen, to_pymatgen

from ..build import MaterialWithBuildMetadata
from ..third_party import PymatgenSpacegroupAnalyzer


class LatticeMaterialAnalyzer(BaseMaterialAnalyzer):
    @property
    def spacegroup_analyzer(self):
        return PymatgenSpacegroupAnalyzer(to_pymatgen(self.material))

    @property
    def material_with_primitive_lattice(self: MaterialWithBuildMetadata) -> MaterialWithBuildMetadata:
        """
        Convert a structure to its primitive cell while preserving atom labels.
        """
        primitive_material = MaterialWithBuildMetadata.create(
            from_pymatgen(self.spacegroup_analyzer.get_primitive_standard_structure())
        )

        # Preserve labels from original material
        return self._preserve_labels_after_lattice_transformation(
            original_material=self.material, transformed_material=primitive_material
        )

    @property
    def material_with_conventional_lattice(self: MaterialWithBuildMetadata) -> MaterialWithBuildMetadata:
        """
        Convert a structure to its conventional cell while preserving atom labels.
        """
        conventional_material = MaterialWithBuildMetadata.create(
            from_pymatgen(self.spacegroup_analyzer.get_conventional_standard_structure())
        )

        # Preserve labels from original material
        return self._preserve_labels_after_lattice_transformation(
            original_material=self.material, transformed_material=conventional_material
        )

    def _preserve_labels_after_lattice_transformation(
        self, original_material: MaterialWithBuildMetadata, transformed_material: MaterialWithBuildMetadata
    ) -> MaterialWithBuildMetadata:
        """
        Preserve atom labels after lattice transformation by mapping atoms based on coordinates.

        Maps each atom in the transformed material to its closest counterpart in the original material
        based on cartesian coordinates and transfers their labels.

        Args:
            original_material: The original material with labels
            transformed_material: The transformed material without labels

        Returns:
            MaterialWithBuildMetadata: Transformed material with restored labels
        """
        if not original_material.basis.labels.values:
            return transformed_material

        labeled_transformed = transformed_material.clone()

        original_cartesian = original_material.clone()
        original_cartesian.to_cartesian()
        original_coords = np.array(original_cartesian.coordinates_array)
        original_labels = original_material.basis.labels.values

        # Get transformed material coordinates in cartesian
        transformed_cartesian = labeled_transformed.clone()
        transformed_cartesian.to_cartesian()
        transformed_coords = np.array(transformed_cartesian.coordinates_array)

        # Map each transformed atom to closest original atom and transfer label
        transformed_labels = []
        for transformed_coord in transformed_coords:
            # Find closest atom in original material
            distances = np.linalg.norm(original_coords - transformed_coord, axis=1)
            closest_original_idx = int(np.argmin(distances))

            # Transfer the label from closest original atom
            closest_label = (
                original_labels[closest_original_idx]
                if closest_original_idx < len(original_labels)
                else original_labels[0]
                if original_labels
                else 0
            )
            transformed_labels.append(closest_label)

        labeled_transformed.set_labels_from_list(transformed_labels)

        return labeled_transformed


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
