from typing import Dict, Optional

import numpy as np
from mat3ra.made.material import Material
from .enums import NanoparticleShapes
from ...build import BaseConfiguration


class NanoparticleConfiguration(BaseConfiguration):
    """
    Configuration for building a nanoparticle.

    Attributes:
        material (Material): The base material for the nanoparticle.
        shape (NanoparticleShapes): The desired shape of the nanoparticle.
        parameters (dict): Dictionary of parameters to pass to the corresponding ASE constructor.
        vacuum_padding (float): Vacuum padding around the nanoparticle.
    """

    material: Material
    shape: NanoparticleShapes
    parameters: Optional[Dict] = None  # Shape-specific parameters (e.g., layers, size)
    vacuum_padding: float = 10.0  # Padding for vacuum space around the nanoparticle

    @property
    def lattice_type(self) -> str:
        return self.material.lattice.type

    @property
    def lattice_constant(self) -> float:
        lattice_constants = [self.material.lattice.a, self.material.lattice.b, self.material.lattice.c]
        # If lattice constants are not equal within a tolerance, raise an error
        if not np.all(np.isclose(lattice_constants, lattice_constants[0], atol=1e-6)):
            raise ValueError("Lattice constants must be equal for isotropic materials")
        lattice_constant = lattice_constants[0]
        return lattice_constant

    @property
    def element(self) -> str:
        return self.material.basis.elements[0]

    @property
    def _json(self):
        return {
            "material": self.material.to_json(),
            "shape": self.shape.value,
            "parameters": self.parameters,
            "vacuum_padding": self.vacuum_padding,
        }
