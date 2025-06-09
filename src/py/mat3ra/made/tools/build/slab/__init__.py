from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from .builders import (
    SlabBuilder,
    AtomicLayersUniqueRepeatedBuilder,
    CrystalLatticePlanesBuilder,
    SlabBuilderParameters,
)
from .configuration import (
    SlabConfiguration,
    CrystalLatticePlanesConfiguration,
    AtomicLayersUniqueRepeatedConfiguration,
    VacuumConfiguration,
)
from ...analyze.termination import Termination
from ...analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer

DEFAULT_XY_SUPERCELL_MATRIX = ([1, 0], [0, 1])


def create_slab(
    crystal: Material,
    miller_indices: tuple[int, int, int] = (0, 0, 1),
    use_conventional_cell=True,
    use_orthogonal_c: bool = True,
    termination: Termination = None,
    number_of_layers=1,
    vacuum=10.0,
    xy_supercell_matrix=DEFAULT_XY_SUPERCELL_MATRIX,
) -> Material:

    crystal_lattice_planes_configuration = CrystalLatticePlanesConfiguration(
        crystal=crystal, miller_indices=miller_indices
    )
    crystal_lattice_planes_material = CrystalLatticePlanesBuilder().get_material(crystal_lattice_planes_configuration)
    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(crystal_lattice_planes_material)

    terminations = crystal_lattice_planes_analyzer.terminations

    if termination is None:
        termination = crystal_lattice_planes_analyzer.default_termination
        print(f"No termination provided. Using default termination: {termination}")
    elif termination not in terminations:
        raise ValueError(f"Termination {termination} not found in available terminations: {terminations}")

    atomic_layers_unique_repeated_configuration = AtomicLayersUniqueRepeatedConfiguration(
        crystal=crystal_lattice_planes_configuration.crystal,
        miller_indices=miller_indices,
        termination_top=termination,
        number_of_repetitions=number_of_layers,
    )

    atomic_layers_unique_repeated_material = AtomicLayersUniqueRepeatedBuilder().get_material(
        atomic_layers_unique_repeated_configuration
    )

    vacuum_configuration = VacuumConfiguration(
        size=vacuum, crystal=atomic_layers_unique_repeated_material, direction=AxisEnum.z
    )

    slab_configuration = SlabConfiguration(
        stack_components=[atomic_layers_unique_repeated_configuration, vacuum_configuration],
        direction=AxisEnum.z,
    )

    slab_builder_parameters = SlabBuilderParameters(
        xy_supercell_matrix=xy_supercell_matrix,
        use_orthogonal_c=use_orthogonal_c,
    )
    builder = SlabBuilder(build_parameters=slab_builder_parameters)
    return builder.get_material(slab_configuration)


def create_slab_if_not(material: Material, default_slab_configuration: SlabConfiguration) -> Material:
    slab = material
    if not slab.metadata or slab.metadata["build"]["configuration"]["type"] != SlabConfiguration.__name__:
        print("The material is not a slab. Creating a new slab...")
        slab = create_slab(default_slab_configuration)
    return slab


def get_slab_terminations(material: Material, miller_indices: tuple[int, int, int] = (0, 0, 1)) -> list[Termination]:
    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
        material=material, miller_indices=miller_indices
    )
    return crystal_lattice_planes_analyzer.terminations


def select_slab_termination(terminations: list[Termination], formula: str = None) -> Termination:
    if not terminations:
        raise ValueError("No terminations available.")
    if formula is None:
        return terminations[0]
    for termination in terminations:
        if termination.formula == formula:
            return termination
    raise ValueError(f"Termination with formula {formula} not found in available terminations: {terminations}")


def get_slab_material_in_standard_representation(slab_material: Material) -> Material:
    """
    Get the slab material in a standard representation.

    Args:
        slab_material (Material): The slab material to convert.

    Returns:
        Material: The slab material in standard representation.
    """
    pass
