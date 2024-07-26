from typing import List, Optional, Any

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder

from .configuration import DeformationConfiguration


class DeformationBuilder(BaseBuilder):
    _ConfigurationType: type(DeformationConfiguration) = DeformationConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    @staticmethod
    def deform_slab(configuration):
        new_material = configuration.slab.clone()
        new_material.to_cartesian()
        new_coordinates = []
        for coord in new_material.basis.coordinates.values:
            perturbed_coord = configuration.deformation_function[0](coord)
            new_coordinates.append(perturbed_coord)
        new_basis = new_material.basis.copy()
        new_basis.coordinates.values = new_coordinates
        new_basis.to_crystal()
        new_material.basis = new_basis
        return new_material

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        """Generate materials with applied deformation based on the given configuration."""
        new_material = self.deform_slab(configuration)
        return [new_material]

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        deformation_details = f"Deformation: {configuration.deformation_function[0].__name__}"
        material.name = f"{material.name} ({deformation_details})"
        return material


class ContinuousDeformationBuilder(DeformationBuilder):
    pass

