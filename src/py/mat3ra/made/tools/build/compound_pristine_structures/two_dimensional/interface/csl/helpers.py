from typing import List, Optional, Tuple, Union

from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from .. import InterfaceBuilder, InterfaceConfiguration
from ......build_components.entities.core.two_dimensional.vacuum.configuration import VacuumConfiguration
from .....pristine_structures.two_dimensional.slab_strained_supercell.configuration import (
    SlabStrainedSupercellConfiguration,
)
from ......analyze.interface import CSLInterfaceAnalyzer
from ......analyze.lattice import get_material_with_conventional_lattice
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.configuration import SlabConfiguration
from ......build_components.metadata import MaterialWithBuildMetadata


def get_csl_strained_configurations(
    substrate_material: Union[Material, MaterialWithBuildMetadata],
    film_material: Union[Material, MaterialWithBuildMetadata],
    substrate_miller_indices: Tuple[int, int, int],
    film_miller_indices: Tuple[int, int, int],
    number_of_layers: int,
    vacuum: float,
    max_area: float = 100.0,
    length_tolerance: float = 0.05,
    angle_step: float = 5.0,
    max_rotation_angle: float = 180.0,
    max_supercell_size: int = 10,
    strain_tolerance: float = 10.0,
    match_id: int = 0,
    use_conventional_cell: bool = True,
) -> Tuple[List[SlabStrainedSupercellConfiguration], float]:
    """
    Get strained configurations for CSL (Coincidence Site Lattice) interface matching.

    This function creates strained configurations for both substrate and film phases
    based on CSL matching with rotational supercells and diagonal substrate matching.

    Args:
        substrate_material (Material): The substrate material.
        film_material (Material): The film material.
        substrate_miller_indices (Tuple[int, int, int]): Miller indices for substrate slab surface.
        film_miller_indices (Tuple[int, int, int]): Miller indices for film slab surface.
        number_of_layers (int): Number of atomic layers in the slabs.
        vacuum (float): Size of the vacuum layer in Angstroms.
        max_area (float): Maximum area for supercell matching.
        length_tolerance (float): Tolerance for matching lattice vector lengths.
        angle_step (float): Step size for rotation angles in degrees.
        max_rotation_angle (float): Maximum rotation angle to test in degrees.
        max_supercell_size (int): Maximum supercell matrix element size.
        strain_tolerance (float): Maximum acceptable strain percentage.
        match_id (int): ID of the match to use (0 for first match).
        use_conventional_cell (bool): Whether to use conventional cells.

    Returns:
        Tuple[List[SlabStrainedSupercellConfiguration], float]:
            List of strained configurations [substrate, film] and the strain percentage.

    Raises:
        ValueError: If no CSL matches are found.
    """
    if use_conventional_cell:
        substrate_material = get_material_with_conventional_lattice(substrate_material)
        film_material = get_material_with_conventional_lattice(film_material)

    substrate_slab_config = SlabConfiguration.from_parameters(
        material_or_dict=substrate_material,
        miller_indices=substrate_miller_indices,
        number_of_layers=number_of_layers,
        vacuum=vacuum,
    )

    film_slab_config = SlabConfiguration.from_parameters(
        material_or_dict=film_material,
        miller_indices=film_miller_indices,
        number_of_layers=number_of_layers,
        vacuum=vacuum,
    )

    analyzer = CSLInterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
        max_area=max_area,
        length_tolerance=length_tolerance,
        angle_step=angle_step,
        max_rotation_angle=max_rotation_angle,
        max_supercell_size=max_supercell_size,
        strain_tolerance=strain_tolerance,
    )

    match_holders = analyzer.csl_match_holders
    if not match_holders:
        raise ValueError("No CSL matches found with the given parameters")

    selected_config = analyzer.get_strained_configuration_by_match_id(match_id)
    strain_percentage = match_holders[match_id].total_strain_percentage
    return [selected_config.substrate_configuration, selected_config.film_configuration], strain_percentage


def create_interface_csl(
    substrate_material: Union[Material, MaterialWithBuildMetadata],
    film_material: Union[Material, MaterialWithBuildMetadata],
    substrate_miller_indices: Tuple[int, int, int] = (0, 0, 1),
    film_miller_indices: Tuple[int, int, int] = (0, 0, 1),
    number_of_layers: int = 1,
    vacuum: float = 0.0,
    direction: AxisEnum = AxisEnum.z,
    gap: float = 3.0,
    max_area: float = 100.0,
    length_tolerance: float = 0.05,
    angle_step: float = 5.0,
    max_rotation_angle: float = 180.0,
    max_supercell_size: int = 10,
    strain_tolerance: float = 10.0,
    match_id: int = 0,
    use_conventional_cell: bool = True,
    remove_overlapping_atoms: bool = True,
    tolerance_for_overlap: float = 1.0,
) -> MaterialWithBuildMetadata:
    """
    Create a CSL (Coincidence Site Lattice) interface between two materials.

    This function creates a CSL interface by:
    1. Creating slab configurations for both materials
    2. Finding CSL matches using rotational supercells and diagonal substrate matching
    3. Creating strained configurations with directional gaps for both phases
    4. Stacking them along the specified direction

    Args:
        substrate_material (Material): The substrate material.
        film_material (Material): The film material.
        substrate_miller_indices (Tuple[int, int, int]): Miller indices for substrate slab surface.
        film_miller_indices (Tuple[int, int, int]): Miller indices for film slab surface.
        number_of_layers (int): Number of atomic layers in the slabs.
        vacuum (float): Size of the vacuum layer in Angstroms.
        direction (AxisEnum): Direction along which to stack components (x, y, or z).
        gap (float): The gap between the two phases in Angstroms.
        max_area (float): Maximum area for supercell matching.
        length_tolerance (float): Tolerance for matching lattice vector lengths.
        angle_step (float): Step size for rotation angles in degrees.
        max_rotation_angle (float): Maximum rotation angle to test in degrees.
        max_supercell_size (int): Maximum supercell matrix element size.
        strain_tolerance (float): Maximum acceptable strain percentage.
        match_id (int): ID of the match to use (0 for first match).
        use_conventional_cell (bool): Whether to use conventional cells.
        remove_overlapping_atoms (bool): Whether to resolve overlapping atoms in the interface after creation.
        tolerance_for_overlap (float): Tolerance for resolving overlapping atoms, in Angstroms.

    Returns:
        Material: The CSL interface material.

    Raises:
        ValueError: If no CSL matches are found.
    """
    strained_configs, _ = get_csl_strained_configurations(
        substrate_material=substrate_material,
        film_material=film_material,
        substrate_miller_indices=substrate_miller_indices,
        film_miller_indices=film_miller_indices,
        number_of_layers=number_of_layers,
        vacuum=vacuum,
        max_area=max_area,
        length_tolerance=length_tolerance,
        angle_step=angle_step,
        max_rotation_angle=max_rotation_angle,
        max_supercell_size=max_supercell_size,
        strain_tolerance=strain_tolerance,
        match_id=match_id,
        use_conventional_cell=use_conventional_cell,
    )

    vacuum_config = VacuumConfiguration(size=vacuum, crystal=None, direction=direction)
    interface_config = InterfaceConfiguration(
        stack_components=strained_configs + [vacuum_config],
        gaps=ArrayWithIds.from_values([gap, gap]),
        direction=direction,
    )

    builder = InterfaceBuilder()
    interface = builder.get_material(interface_config)
    if remove_overlapping_atoms:
        interface.basis.resolve_colliding_coordinates(tolerance=tolerance_for_overlap)
    return interface
