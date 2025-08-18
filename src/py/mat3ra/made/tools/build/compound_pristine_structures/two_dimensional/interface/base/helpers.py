from typing import List, Optional, Tuple, Union

from mat3ra.code.array_with_ids import ArrayWithIds

from mat3ra.made.material import Material
from .build_parameters import InterfaceBuilderParameters
from .builder import InterfaceBuilder
from .configuration import InterfaceConfiguration
from ......analyze.interface import InterfaceAnalyzer
from ......analyze.slab import SlabMaterialAnalyzer
from ......build.pristine_structures.two_dimensional.slab.helpers import create_slab
from ......build_components import MaterialWithBuildMetadata
from ......build_components.entities.core.two_dimensional.vacuum.configuration import VacuumConfiguration


def create_interface_simple_between_slabs(
    substrate_slab: MaterialWithBuildMetadata,
    film_slab: MaterialWithBuildMetadata,
    gap: Optional[float] = None,
    vacuum: float = 10.0,
    xy_shift: List[float] = [0, 0],
    reduce_result_cell_to_primitive: bool = False,
) -> Material:
    """
    Create an interface between two slab materials with specified parameters.
    No strain matching is performed, the interface is created as a simple stack of two slabs.

    Args:
        substrate_slab (MaterialWithBuildMetadata): Substrate slab material.
        film_slab (MaterialWithBuildMetadata): Film slab material.
        gap (Optional[float]): Gap between the two materials, in Angstroms.
            Distance between top most atom of the substrate and bottom most atom of the film.
        vacuum (float): Size of the vacuum layer in Angstroms.
        xy_shift (List[float]): Shift in x and y directions, in Angstroms.

    Returns:
        Material: The interface material.
    """
    substrate_analyzer = SlabMaterialAnalyzer(material=substrate_slab)
    film_analyzer = SlabMaterialAnalyzer(material=film_slab)

    analyzer = InterfaceAnalyzer(
        substrate_slab_configuration=substrate_analyzer.build_configuration,
        film_slab_configuration=film_analyzer.build_configuration,
        substrate_build_parameters=substrate_analyzer.build_parameters,
        film_build_parameters=film_analyzer.build_parameters,
    )

    film_configuration = analyzer.film_strained_configuration
    substrate_configuration = analyzer.substrate_strained_configuration
    vacuum_configuration = VacuumConfiguration(
        size=vacuum,
    )

    config = InterfaceConfiguration(
        stack_components=[substrate_configuration, film_configuration, vacuum_configuration],
        xy_shift=xy_shift,
        gaps=ArrayWithIds.from_values([gap]),
    )
    builder = InterfaceBuilder(
        build_parameters=InterfaceBuilderParameters(make_primitive=reduce_result_cell_to_primitive)
    )
    interface = builder.get_material(config)
    return interface


def create_interface_simple(
    substrate_crystal: Union[Material, MaterialWithBuildMetadata],
    film_crystal: Union[Material, MaterialWithBuildMetadata],
    substrate_miller_indices: Tuple[int, int, int] = (0, 0, 1),
    film_miller_indices: Tuple[int, int, int] = (0, 0, 1),
    substrate_number_of_layers: int = 1,
    film_number_of_layers: int = 1,
    substrate_termination_formula: Optional[str] = None,
    film_termination_formula: Optional[str] = None,
    substrate_xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]],
    film_xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]],
    gap: Optional[float] = None,
    vacuum: float = 10.0,
    xy_shift: List[float] = [0, 0],
    use_conventional_cell: bool = True,
    reduce_result_cell_to_primitive: bool = False,
) -> Material:
    """
    Create an interface between two crystal materials by first creating slabs and then stacking them.
    No strain matching is performed, the interface is created as a simple stack of two slabs.

    Args:
        substrate_crystal (Union[Material, MaterialWithBuildMetadata]): Substrate crystal material.
        film_crystal (Union[Material, MaterialWithBuildMetadata]): Film crystal material.
        substrate_miller_indices (Tuple[int, int, int]): Miller indices for the substrate slab surface.
        film_miller_indices (Tuple[int, int, int]): Miller indices for the film slab surface.
        substrate_number_of_layers (int): Number of atomic layers in the substrate slab.
        film_number_of_layers (int): Number of atomic layers in the film slab.
        substrate_termination_formula (Optional[str]): Formula of the termination for substrate slab.
        film_termination_formula (Optional[str]): Formula of the termination for film slab.
        substrate_xy_supercell_matrix (List[List[int]]): Supercell matrix for substrate slab xy plane.
        film_xy_supercell_matrix (List[List[int]]): Supercell matrix for film slab xy plane.
        gap (Optional[float]): Gap between the two materials, in Angstroms.
        vacuum (float): Size of the vacuum layer in Angstroms.
        xy_shift (List[float]): Shift in x and y directions, in Angstroms.
        use_conventional_cell (bool): Whether to use conventional cell for crystals.
        reduce_result_cell_to_primitive (bool): Whether to reduce result cell to primitive.

    Returns:
        Material: The interface material.
    """
    substrate_slab = create_slab(
        crystal=substrate_crystal,
        miller_indices=substrate_miller_indices,
        number_of_layers=substrate_number_of_layers,
        vacuum=0,
        termination_top_formula=substrate_termination_formula,
        xy_supercell_matrix=substrate_xy_supercell_matrix,
        use_conventional_cell=use_conventional_cell,
    )
    film_slab = create_slab(
        crystal=film_crystal,
        miller_indices=film_miller_indices,
        number_of_layers=film_number_of_layers,
        vacuum=0,
        termination_top_formula=film_termination_formula,
        xy_supercell_matrix=film_xy_supercell_matrix,
        use_conventional_cell=use_conventional_cell,
    )

    return create_interface_simple_between_slabs(
        substrate_slab=substrate_slab,
        film_slab=film_slab,
        gap=gap,
        vacuum=vacuum,
        xy_shift=xy_shift,
        reduce_result_cell_to_primitive=reduce_result_cell_to_primitive,
    )
