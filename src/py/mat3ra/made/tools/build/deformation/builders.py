from typing import List, Optional, Any

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder

from scipy.integrate import quad
from scipy.optimize import root_scalar
from .configuration import DeformationConfiguration


class DeformationBuilder(BaseBuilder):
    _ConfigurationType: type(DeformationConfiguration) = DeformationConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    @staticmethod
    def deform_slab(configuration):
        new_material = configuration.slab.clone()
        # new_material.to_cartesian()
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
    def df_dx(self, x, wavelength, phase):
        return (2 * np.pi / wavelength) * np.cos(2 * np.pi * x / wavelength + phase)

    def arc_length_integral(self, x_prime, x, amplitude, wavelength, phase):
        def integrand(t):
            return np.sqrt(1 + (amplitude * 2 * np.pi / wavelength * np.cos(2 * np.pi * t / wavelength + phase)) ** 2)

        # Compute the arc length from 0 to x_prime.
        calculated_length = quad(func=integrand, a=0, b=x_prime)[0]
        return calculated_length - x

    def find_x_prime(self, x, amplitude, wavelength, phase):
        # Find x' such that the integral from 0 to x' equals x
        result = root_scalar(
            self.arc_length_integral, args=(x, amplitude, wavelength, phase), bracket=[0, 10 * x], method="brentq"
        )
        return result.root

    def deform_slab_continuously(self, configuration):
        new_material = configuration.slab.clone()
        new_material.to_cartesian()
        new_coordinates = []
        for coord in new_material.basis.coordinates.values:
            x_prime = self.find_x_prime(
                coord[0],
                configuration.deformation_function[1]["amplitude"],
                configuration.deformation_function[1]["wavelength"],
                configuration.deformation_function[1]["phase"],
            )
            print(coord[0], x_prime)
            perturbed_coord = configuration.deformation_function[0]([x_prime, coord[1], coord[2]])
            new_coordinates.append(perturbed_coord)

        new_basis = new_material.basis.copy()
        new_basis.coordinates.values = new_coordinates
        new_basis.to_crystal()
        new_material.basis = new_basis

        return new_material

    def _generate(
        self, configuration: DeformationBuilder._ConfigurationType
    ) -> List[DeformationBuilder._GeneratedItemType]:
        """Generate materials with applied continuous deformation based on the given configuration."""
        new_material = self.deform_slab_continuously(configuration)
        return [new_material]
