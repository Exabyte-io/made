from typing import Dict, Optional, Tuple

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.tools.convert import decorator_convert_material_args_kwargs_to_structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ...convert import from_pymatgen
from ...third_party import PymatgenStructure

from .enums import NanoparticleShapes
from ...build import BaseConfiguration


class BaseNanoparticleConfiguration(BaseConfiguration):
    material: Material  # Material to use for the nanoparticle
    vacuum_padding: float = 10.0  # Padding for vacuum space around the nanoparticle

    @property
    def lattice_type(self) -> str:
        return self.material.lattice.type

    @property
    def lattice_constant(self) -> float:
        material = self.material
        conventional_material = ASEBasedNanoparticleConfiguration.convert_to_conventional(material)
        lattice_constants = [
            conventional_material.lattice.a,
            conventional_material.lattice.b,
            conventional_material.lattice.c,
        ]
        # If lattice constants are not equal within a tolerance, raise an error
        if not np.all(np.isclose(lattice_constants, lattice_constants[0], atol=1e-6)):
            raise ValueError("Lattice constants must be equal for isotropic materials")
        lattice_constant = lattice_constants[0]
        return lattice_constant

    @property
    def element(self) -> str:
        return self.material.basis.elements.get_element_value_by_index(0)

    @property
    def _json(self):
        return {
            "material": self.material.to_json(),
            "vacuum_padding": self.vacuum_padding,
        }

    # TODO: move to a separate module
    @staticmethod
    @decorator_convert_material_args_kwargs_to_structure
    def convert_to_primitive(structure: PymatgenStructure) -> Material:
        """
        Convert a structure to its primitive cell.
        """
        analyzer = SpacegroupAnalyzer(structure)
        return Material(from_pymatgen(analyzer.get_primitive_standard_structure()))

    @staticmethod
    @decorator_convert_material_args_kwargs_to_structure
    def convert_to_conventional(structure: PymatgenStructure) -> Material:
        """
        Convert a structure to its conventional cell.
        """
        analyzer = SpacegroupAnalyzer(structure)
        return Material(from_pymatgen(analyzer.get_conventional_standard_structure()))


class NanoparticleConfiguration(BaseNanoparticleConfiguration):
    shape: NanoparticleShapes = NanoparticleShapes.ICOSAHEDRON  # Shape of the nanoparticle
    orientation_z: Tuple[int, int, int] = (0, 0, 1)  # Orientation of the crystallographic axis in the z-direction
    radius: float = 5.0  # Radius of the nanoparticle (largest feature size for a shape), in Angstroms

    @property
    def _json(self):
        return {
            **super()._json,
            "shape": self.shape.value,
            "orientation_z": self.orientation_z,
            "radius": self.radius,
        }


class ASEBasedNanoparticleConfiguration(BaseNanoparticleConfiguration):
    """
    Configuration for building a nanoparticle.

    Attributes:
        material (Material): The base material for the nanoparticle.
        shape (NanoparticleShapes): The desired shape of the nanoparticle.
        parameters (dict): Dictionary of parameters to pass to the corresponding ASE constructor.
        vacuum_padding (float): Vacuum padding around the nanoparticle.
    """

    shape: NanoparticleShapes
    parameters: Optional[Dict] = None  # Shape-specific parameters (e.g., layers, size)

    @property
    def _json(self):
        return {
            **super()._json,
            "shape": self.shape.value,
            "parameters": self.parameters,
        }
