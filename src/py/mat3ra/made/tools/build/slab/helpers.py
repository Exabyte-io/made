from typing import List, Tuple, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.entities import Termination
from .builders import (
    SlabBuilder,
    SlabBuilderParameters,
    AtomicLayersUniqueRepeatedBuilder,
)
from .configurations import (
    SlabConfiguration,
    AtomicLayersUniqueRepeatedConfiguration,
)
from ..metadata import MaterialWithBuildMetadata
from ...analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer

DEFAULT_XY_SUPERCELL_MATRIX = ([1, 0], [0, 1])


def create_atomic_layers(
    material: Union[Material, MaterialWithBuildMetadata],
    miller_indices: Tuple[int, int, int] = (0, 0, 1),
    termination: Termination = None,
    number_of_layers: int = 1,
) -> Material:
    """
    Creates a material composed of repeated unique atomic layers from a given crystal.

    This function identifies the sequence of unique atomic layers for the given Miller
    indices and then constructs the material by repeating this sequence, starting
    with a given surface termination.

    Args:
        material (Material): The crystal material to create atomic layers from.
        miller_indices (Tuple[int, int, int]): Miller indices for the atomic layers.
        termination (Termination): The termination to use for the atomic layers.
        number_of_layers (int): Number of times to repeat the sequence of unique atomic layers.
    Returns:
        Material: The atomic layers material.

    """
    atomic_layers_config = AtomicLayersUniqueRepeatedConfiguration(
        crystal=material,
        miller_indices=miller_indices,
        termination_top=termination,
        number_of_repetitions=number_of_layers,
    )

    atomic_layers_builder = AtomicLayersUniqueRepeatedBuilder()
    atomic_layers_material = atomic_layers_builder.get_material(atomic_layers_config)

    return atomic_layers_material


def create_slab(
    crystal: Union[Material, MaterialWithBuildMetadata],
    miller_indices: Tuple[int, int, int] = (0, 0, 1),
    use_conventional_cell=True,
    use_orthogonal_c: bool = True,
    termination: Termination = None,
    termination_formula: str = None,
    number_of_layers=1,
    vacuum=10.0,
    xy_supercell_matrix=DEFAULT_XY_SUPERCELL_MATRIX,
) -> MaterialWithBuildMetadata:
    """
    Creates a slab material from a crystal material with specified Miller indices and other parameters.

    Args:
        crystal (Material): The crystal material to create the slab from.
        miller_indices (Tuple[int, int, int]): Miller indices for the slab surface.
        use_conventional_cell (bool): Whether to use the conventional cell for the crystal to apply Miller indices.
        use_orthogonal_c (bool): Whether to make slab with c-lattice orthogonal to miller indices plane, along z-axis.
        termination (Termination): The termination to use for the slab, supply a class instance of Termination.
        termination_formula (str): The stoichiometric formula of the termination, e.g. "Si", "TiO", "SrTiO", "Hf2O".
        number_of_layers (int): Number of atomic layers in the slab, in the number of unit cells.
        vacuum (float): Size of the vacuum layer in Angstroms.
        xy_supercell_matrix (List[List[int]]): Supercell matrix for the xy plane to apply to the generated slab.

    Returns:
        Material: The generated slab material.
    """
    material_to_use = crystal

    if use_conventional_cell:
        crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
            material=crystal, miller_indices=miller_indices
        )
        material_to_use = crystal_lattice_planes_analyzer.material_with_conventional_lattice

    if termination is not None:
        termination_formula = termination.formula

    slab_builder_parameters = SlabBuilderParameters(
        xy_supercell_matrix=xy_supercell_matrix,
        use_orthogonal_c=use_orthogonal_c,
    )
    slab_configuration = SlabConfiguration.from_parameters(
        material_or_dict=material_to_use,
        miller_indices=miller_indices,
        number_of_layers=number_of_layers,
        termination_formula=termination_formula,
        vacuum=vacuum,
        use_conventional_cell=use_conventional_cell,
    )
    builder = SlabBuilder(build_parameters=slab_builder_parameters)
    return builder.get_material(slab_configuration)


def create_slab_if_not(
    material: MaterialWithBuildMetadata, default_slab_configuration: SlabConfiguration
) -> MaterialWithBuildMetadata:
    slab = material
    build_metadata_slab = slab.metadata.get_build_metadata_of_type("SlabConfiguration")
    if build_metadata_slab is None:
        print("The material is not a slab. Creating a new slab...")
        slab = SlabBuilder().get_material(default_slab_configuration)
    return slab


def get_slab_terminations(
    material: Union[Material, MaterialWithBuildMetadata], miller_indices: Tuple[int, int, int] = (0, 0, 1)
) -> List[Termination]:
    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
        material=material, miller_indices=miller_indices
    )
    return crystal_lattice_planes_analyzer.terminations


def get_slab_material_in_standard_representation(slab_material: Union[Material, MaterialWithBuildMetadata]) -> Material:
    """
    Get the slab material in a standard representation.

    Args:
        slab_material (Material): The slab material to convert.

    Returns:
        Material: The slab material in standard representation.
    """
    pass
