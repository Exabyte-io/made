from typing import Callable, List, Optional, Union

import numpy as np
from mat3ra.made.material import Material
from pydantic import BaseModel

from ..analyze import get_surface_atom_indices
from ..convert import from_ase
from ..convert.utils import InterfacePartsEnum
from ..enums import SurfaceTypes
from ..modify import get_interface_part
from ..third_party import ASEAtoms, ASECalculator, ASEFixAtoms, ASEFixedPlane, ase_all_changes
from .interaction_functions import sum_of_inverse_distances_squared


class MaterialCalculatorParameters(BaseModel):
    """
    Defines the parameters for a material calculator.

    Args:
        interaction_function (Callable): A function used to calculate the interaction metric between
        sets of coordinates. The default function is sum_of_inverse_distances_squared.
    """

    interaction_function: Callable = sum_of_inverse_distances_squared


class InterfaceMaterialCalculatorParameters(MaterialCalculatorParameters):
    """
    Parameters specific to the calculation of interaction energies between
    an interface material's film and substrate.

    Args:
        shadowing_radius (float): Radius used to determine the surface atoms of the film or substrate
                                  for interaction calculations. Default is 2.5 Å.
    """

    shadowing_radius: float = 2.5


class MaterialCalculator(BaseModel):
    """
    A base class for performing calculations on materials.

    This class uses the parameters defined in MaterialCalculatorParameters to calculate
    interaction metrics between atoms or sets of coordinates within the material.

    Args:
        calculator_parameters (MaterialCalculatorParameters): Parameters controlling the calculator,
        including the interaction function.
    """

    calculator_parameters: MaterialCalculatorParameters = MaterialCalculatorParameters()

    def get_energy(self, material: Material):
        """
        Calculate the energy (or other metric) for a material.

        Args:
            material (Material): The material to calculate the interaction energy for.

        Returns:
            float: The interaction energy between the coordinates of the material,
            calculated using the specified interaction function.
        """
        return self.calculator_parameters.interaction_function(material.coordinates, material.coordinates)


class InterfaceMaterialCalculator(MaterialCalculator):
    """
    A specialized calculator for computing the interaction energy between a film and a substrate
    in an interface material.

    This class extends MaterialCalculator and uses additional parameters specific to interface materials,
    such as the shadowing radius to detect surface atoms.

    Args:
        calculator_parameters (InterfaceMaterialCalculatorParameters): Parameters that include the
        shadowing radius and interaction function.
    """

    calculator_parameters: InterfaceMaterialCalculatorParameters = InterfaceMaterialCalculatorParameters()

    def get_energy(
        self,
        material: Material,
        shadowing_radius: float = 2.5,
        interaction_function: Callable = sum_of_inverse_distances_squared,
    ) -> float:
        """
        Calculate the interaction energy between the film and substrate in an interface material.

        This method uses the shadowing radius to detect surface atoms and applies the given
        interaction function to calculate the interaction between the film and substrate.

        Args:
            material (Material): The interface Material object consisting of both the film and substrate.
            shadowing_radius (float): The radius used to detect surface atoms for the interaction
                                      calculation. Defaults to 2.5 Å.
            interaction_function (Callable): A function to compute the interaction between the film and
                                             substrate. Defaults to sum_of_inverse_distances_squared.

        Returns:
            float: The calculated interaction energy between the film and substrate.
        """
        film_material = get_interface_part(material, part=InterfacePartsEnum.FILM)
        substrate_material = get_interface_part(material, part=InterfacePartsEnum.SUBSTRATE)

        film_surface_atom_indices = get_surface_atom_indices(
            film_material, SurfaceTypes.BOTTOM, shadowing_radius=shadowing_radius
        )
        substrate_surface_atom_indices = get_surface_atom_indices(
            substrate_material, SurfaceTypes.TOP, shadowing_radius=shadowing_radius
        )

        film_surface_atom_coordinates = film_material.basis.coordinates
        film_surface_atom_coordinates.filter_by_ids(film_surface_atom_indices)
        substrate_surface_atom_coordinates = substrate_material.basis.coordinates
        substrate_surface_atom_coordinates.filter_by_ids(substrate_surface_atom_indices)

        film_coordinates_values = np.array(film_surface_atom_coordinates.values)
        substrate_coordinates_values = np.array(substrate_surface_atom_coordinates.values)

        return interaction_function(film_coordinates_values, substrate_coordinates_values)


class SurfaceDistanceCalculator(ASECalculator):
    """
    ASE calculator that computes the norm of distances between interfacial gap facing atoms
    of the film and the substrate.

    Args:
        shadowing_radius (float): The radius for atom to shadow underlying from being considered surface, in Angstroms.
        force_constant (float): The force constant for the finite difference approximation of the forces.
        fix_substrate (bool): Whether to fix the substrate atoms.
        fix_z (bool): Whether to fix atoms movement in the z direction.
        symprec (float): The symmetry precision for the ASE calculator.

    Example usage:
    ```python
    from ase.optimize import BFGS
    atoms = to_ase(material)
    calc = SurfaceDistanceCalculator(shadowing_radius=2.5)

    atoms.calc = calc
    opt = BFGS(atoms)
    opt.run(fmax=0.05)
    ```
    Args:
        shadowing_radius (float): Radius for atoms shadowing underlying from being treated as a surface, in Angstroms.
        force_constant (float): The force constant for the finite difference approximation of the
    Note:
        Built following: https://wiki.fysik.dtu.dk/ase/development/calculators.html

        The calculate method is responsible for computing the energy and forces (if requested).
        Forces are estimated using a finite difference method, which is a simple approximation
        and might not be the most accurate or efficient for all cases.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        shadowing_radius: float = 2.5,
        force_constant: float = 1.0,
        fix_substrate: bool = True,
        fix_z: bool = True,
        symprec: float = 0.01,
        calculator_parameters: MaterialCalculatorParameters = MaterialCalculatorParameters(),
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.shadowing_radius = shadowing_radius
        self.force_constant = force_constant
        self.fix_substrate = fix_substrate
        self.fix_z = fix_z
        self.symprec = symprec
        self.material_calculator = MaterialCalculator(calculator_parameters=calculator_parameters)

    def _add_constraints(self, atoms: ASEAtoms) -> ASEAtoms:
        constraints: List[Union[ASEFixAtoms, ASEFixedPlane]] = []
        if self.fix_substrate:
            substrate_indices = [i for i, tag in enumerate(atoms.get_tags()) if tag == 0]
            constraints.append(ASEFixAtoms(indices=substrate_indices))
        if self.fix_z:
            all_indices = list(range(len(atoms)))
            constraints.append(ASEFixedPlane(indices=all_indices, direction=[0, 0, 1]))

        atoms.set_constraint(constraints)
        return atoms

    def _calculate_forces(self, atoms: ASEAtoms, energy: float) -> np.ndarray:
        forces = np.zeros((len(atoms), 3))
        dx = 0.01
        for i in range(len(atoms)):
            for j in range(3):
                atoms_plus = atoms.copy()
                atoms_plus.positions[i, j] += dx
                material_plus = Material(from_ase(atoms_plus))
                energy_plus = self.material_calculator.get_energy(material_plus)

                forces[i, j] = -self.force_constant * (energy_plus - energy) / dx

        return forces

    def calculate(self, atoms: Optional[ASEAtoms] = None, properties=["energy"], system_changes=ase_all_changes):
        if atoms is None:
            atoms = self.atoms.copy()

        atoms = self._add_constraints(atoms)
        constraints = atoms.constraints

        ASECalculator.calculate(self, atoms, properties, system_changes)
        material = Material(from_ase(atoms))
        energy = self.material_calculator.get_energy(material)

        self.results = {"energy": energy}

        if "forces" in properties:
            forces = self._calculate_forces(atoms, energy)
            for constraint in constraints:
                constraint.adjust_forces(atoms, forces)
            self.results["forces"] = forces

    def get_potential_energy(self, atoms=None, force_consistent=False):
        return self.get_property("energy", atoms)

    def get_forces(self, atoms=None):
        return self.get_property("forces", atoms)
