from typing import List, Optional, Any

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder

from .configuration import PerturbationConfiguration


class PerturbationBuilder(BaseBuilder):
    _ConfigurationType: type(PerturbationConfiguration) = PerturbationConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        """Generate materials with applied perturbations based on the given configuration."""
        new_material = configuration.slab.clone()
        new_material.to_cartesian()
        new_coordinates = []
        for coord in new_material.basis.coordinates.values:
            x, y, z = coord
            perturbed_z = z + configuration.amplitude * configuration.perturbation_func(x, y)
            new_coordinates.append([x, y, perturbed_z])
        new_basis = new_material.basis.copy()
        new_basis.coordinates.values = new_coordinates
        new_basis.to_crystal()
        new_material.basis = new_basis
        return [new_material]

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        """Update material name based on the perturbation details."""
        perturbation_details = f"perturbed_amplitude_{configuration.amplitude}"
        material.name = f"{material.name} ({perturbation_details})"
        return material
