from typing import Callable, List, Optional, Union

import numpy as np
from mat3ra.made.tools.convert.utils import InterfacePartsEnum
from pydantic import BaseModel

from ..material import Material
from .analyze import get_surface_area, get_surface_atom_indices
from .build.interface.utils import get_slab
from .convert import decorator_convert_material_args_kwargs_to_atoms, from_ase
from .enums import SurfaceTypes
from .modify import get_interface_part
from .third_party import ASEAtoms, ASECalculator, ASECalculatorEMT, ASEFixAtoms, ASEFixedPlane, ase_all_changes
from .utils import decorator_handle_periodic_boundary_conditions, get_sum_of_inverse_distances_squared


@decorator_convert_material_args_kwargs_to_atoms
def calculate_total_energy(atoms: ASEAtoms, calculator: ASECalculator):
    """
    Set calculator for ASE Atoms and calculate the total energy.

    Args:
        atoms (ASEAtoms): The Atoms object to calculate the energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The total energy of the atoms.
    """
    atoms.set_calculator(calculator)
    return atoms.get_total_energy()


@decorator_convert_material_args_kwargs_to_atoms
def calculate_total_energy_per_atom(atoms: ASEAtoms, calculator: ASECalculator):
    """
    Set calculator for ASE Atoms and calculate the total energy per atom.

    Args:
            atoms (ASEAtoms): The Atoms object to calculate the energy of.
            calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
            float: The energy per atom of the atoms.
    """
    return calculate_total_energy(atoms, calculator) / atoms.get_global_number_of_atoms()


@decorator_convert_material_args_kwargs_to_atoms
def calculate_surface_energy(slab: ASEAtoms, bulk: ASEAtoms, calculator: ASECalculator):
    """
    Calculate the surface energy by subtracting the weighted bulk energy from the slab energy.

    Args:
        slab (ASEAtoms): The slab Atoms object to calculate the surface energy of.
        bulk (ASEAtoms): The bulk Atoms object to calculate the surface energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The surface energy of the slab.
    """
    number_of_atoms = slab.get_global_number_of_atoms()
    area = get_surface_area(slab)
    return (
        calculate_total_energy(slab, calculator) - calculate_total_energy_per_atom(bulk, calculator) * number_of_atoms
    ) / (2 * area)


@decorator_convert_material_args_kwargs_to_atoms
def calculate_adhesion_energy(
    interface: ASEAtoms, substrate_slab: ASEAtoms, film_slab: ASEAtoms, calculator: ASECalculator
):
    """
    Calculate the adhesion energy.
    The adhesion energy is the difference between the energy of the interface and
    the sum of the energies of the substrate and film.
    According to: 10.1088/0953-8984/27/30/305004

    Args:
        interface (ASEAtoms): The interface ASEAtoms object to calculate the adhesion energy of.
        substrate_slab (ASEAtoms): The substrate slab ASEAtoms object to calculate the adhesion energy of.
        film_slab (ASEAtoms): The film slab ASEAtoms object to calculate the adhesion energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The adhesion energy of the interface.
    """
    energy_substrate_slab = calculate_total_energy(substrate_slab, calculator)
    energy_film_slab = calculate_total_energy(film_slab, calculator)
    energy_interface = calculate_total_energy(interface, calculator)
    area = get_surface_area(interface)
    return (energy_substrate_slab + energy_film_slab - energy_interface) / area


def calculate_interfacial_energy(
    interface: Material,
    substrate_slab: Optional[Material] = None,
    substrate_bulk: Optional[Material] = None,
    film_slab: Optional[Material] = None,
    film_bulk: Optional[Material] = None,
    calculator: ASECalculator = ASECalculatorEMT(),
):
    """
    Calculate the interfacial energy.
    The interfacial energy is the sum of the surface energies of the substrate and film minus the adhesion energy.
    According to Dupré's formula

    Args:
        interface (Material): The interface Material object to calculate the interfacial energy of.
        substrate_slab (Material): The substrate slab Material object to calculate the interfacial energy of.
        substrate_bulk (Material): The substrate bulk Material object to calculate the interfacial energy of.
        film_slab (Material): The film slab Material object to calculate the interfacial energy of.
        film_bulk (Material): The film bulk Material object to calculate the interfacial energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation

    Returns:
        float: The interfacial energy of the interface.
    """
    substrate_slab = get_slab(interface, part="substrate") if substrate_slab is None else substrate_slab
    film_slab = get_slab(interface, part="film") if film_slab is None else film_slab

    build_configuration = interface.metadata["build"]["configuration"] if "build" in interface.metadata else {}
    try:
        substrate_bulk = (
            Material(build_configuration["substrate_configuration"]["bulk"])
            if substrate_bulk is None
            else substrate_bulk
        )
        film_bulk = Material(build_configuration["film_configuration"]["bulk"]) if film_bulk is None else film_bulk
    except KeyError:
        raise ValueError("The substrate and film bulk materials must be provided or defined in the interface metadata.")
    surface_energy_substrate = calculate_surface_energy(substrate_slab, substrate_bulk, calculator)
    surface_energy_film = calculate_surface_energy(film_slab, film_bulk, calculator)
    adhesion_energy = calculate_adhesion_energy(interface, substrate_slab, film_slab, calculator)
    return surface_energy_film + surface_energy_substrate - adhesion_energy


class InteractionCalculatorParameters(BaseModel):
    shadowing_radius: float = 2.5
    interaction_function: Callable = get_sum_of_inverse_distances_squared


@decorator_handle_periodic_boundary_conditions(cutoff=0.25)
def calculate_film_substrate_interaction_metric(
    material: Material,
    shadowing_radius: float = 2.5,
    interaction_function: Callable = get_sum_of_inverse_distances_squared,
) -> float:
    """
    Calculate the interaction metric between the film and substrate.
    Args:
        material (Material): The interface Material object.
        shadowing_radius (float): The shadowing radius to detect the surface atoms, in Angstroms.
        interaction_function (Callable): The metric function to use for the calculation of the interaction.

    Returns:
        float: The calculated norm.
    """
    film_material = get_interface_part(material, part=InterfacePartsEnum.FILM)
    substrate_material = get_interface_part(material, part=InterfacePartsEnum.SUBSTRATE)
    film_atoms_surface_indices = get_surface_atom_indices(
        film_material, SurfaceTypes.BOTTOM, shadowing_radius=shadowing_radius
    )
    substrate_atoms_surface_indices = get_surface_atom_indices(
        substrate_material, SurfaceTypes.TOP, shadowing_radius=shadowing_radius
    )

    film_atoms_surface_coordinates = film_material.basis.coordinates
    film_atoms_surface_coordinates.filter_by_ids(film_atoms_surface_indices)
    substrate_atoms_surface_coordinates = substrate_material.basis.coordinates
    substrate_atoms_surface_coordinates.filter_by_ids(substrate_atoms_surface_indices)

    film_coordinates_values = np.array(film_atoms_surface_coordinates.values)
    substrate_coordinates_values = np.array(substrate_atoms_surface_coordinates.values)

    return interaction_function(film_coordinates_values, substrate_coordinates_values)


class SurfaceDistanceCalculator(ASECalculator):
    """
    ASE calculator that computes the norm of distances between interfacial gap facing atoms
    of the film and the substrate.

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
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.shadowing_radius = shadowing_radius
        self.force_constant = force_constant
        self.fix_substrate = fix_substrate
        self.fix_z = fix_z
        self.symprec = symprec

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

    def _calculate_forces(self, atoms: ASEAtoms, norm: float) -> np.ndarray:
        forces = np.zeros((len(atoms), 3))
        dx = 0.01
        for i in range(len(atoms)):
            for j in range(3):
                atoms_plus = atoms.copy()
                atoms_plus.positions[i, j] += dx
                material_plus = Material(from_ase(atoms_plus))
                norm_plus = calculate_film_substrate_interaction_metric(material_plus, self.shadowing_radius)

                forces[i, j] = -self.force_constant * (norm_plus - norm) / dx

        return forces

    @decorator_convert_material_args_kwargs_to_atoms
    def calculate(self, atoms: Optional[ASEAtoms] = None, properties=["energy"], system_changes=ase_all_changes):
        if atoms is None:
            atoms = self.atoms.copy()

        atoms = self._add_constraints(atoms)
        constraints = atoms.constraints

        ASECalculator.calculate(self, atoms, properties, system_changes)
        material = Material(from_ase(atoms))
        norm = calculate_film_substrate_interaction_metric(material, self.shadowing_radius)

        self.results = {"energy": norm}

        if "forces" in properties:
            forces = self._calculate_forces(atoms, norm)
            for constraint in constraints:
                constraint.adjust_forces(atoms, forces)
            self.results["forces"] = forces

    def get_potential_energy(self, atoms=None, force_consistent=False):
        return self.get_property("energy", atoms)

    def get_forces(self, atoms=None):
        return self.get_property("forces", atoms)
