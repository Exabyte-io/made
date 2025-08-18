from typing import List, Tuple, Union, Optional

from mat3ra.made.material import Material
from . import SlabBuilderParameters, SlabConfiguration, SlabBuilder
from .....analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from .....build_components import MaterialWithBuildMetadata
from .....build_components.entities.auxiliary.two_dimensional.termination import Termination


DEFAULT_XY_SUPERCELL_MATRIX = ([1, 0], [0, 1])


def create_slab(
    crystal: Union[Material, MaterialWithBuildMetadata],
    miller_indices: Tuple[int, int, int] = (0, 0, 1),
    use_conventional_cell=True,
    use_orthogonal_c: bool = True,
    termination_top: Optional[Termination] = None,
    termination_bottom: Optional[Termination] = None,
    termination_top_formula: Optional[str] = None,
    termination_bottom_formula: Optional[str] = None,
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

    if termination_top is not None:
        termination_top_formula = termination_top.formula
    if termination_bottom is not None:
        termination_bottom_formula = termination_bottom.formula

    slab_builder_parameters = SlabBuilderParameters(
        xy_supercell_matrix=xy_supercell_matrix,
        use_orthogonal_c=use_orthogonal_c,
    )
    slab_configuration = SlabConfiguration.from_parameters(
        material_or_dict=material_to_use,
        miller_indices=miller_indices,
        number_of_layers=number_of_layers,
        termination_top_formula=termination_top_formula,
        termination_bottom_formula=termination_bottom_formula,
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
