from mat3ra.code.vector import Vector3D
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.operations.core.unary import translate
from .builders import SlabBuilder, ConventionalCellBuilder, AtomicLayersUniqueRepeatedBuilder
from .configuration import (
    SlabConfiguration,
    ConventionalCellConfiguration,
    CrystalLatticePlanesConfiguration,
    AtomicLayersUniqueRepeatedConfiguration,
    VacuumConfiguration,
    get_terminations,
)
from .utils import select_termination
from .termination import Termination


def create_slab(
    crystal: Material,
    miller_indices: tuple[int, int, int] = (0, 0, 1),
    use_conventional_cell=True,
    use_orthogonal_c: bool = True,
    termination: Termination = None,
    number_of_layers=1,
    vacuum=10.0,
    xy_supercell_matrix=None,
) -> Material:
    if xy_supercell_matrix is None:
        xy_supercell_matrix = [[1, 0], [0, 1]]

    conventional_cell_config = ConventionalCellConfiguration(crystal=crystal)
    working_crystal = (
        ConventionalCellBuilder().get_material(conventional_cell_config) if use_conventional_cell else crystal
    )

    crystal_lattice_planes = CrystalLatticePlanesConfiguration(crystal=working_crystal, miller_indices=miller_indices)

    if termination is None:
        terminations = get_terminations(crystal=working_crystal, miller_indices=miller_indices)
        termination = terminations[0] if terminations else None

    atomic_layers_repeated_config = AtomicLayersUniqueRepeatedConfiguration(
        crystal=crystal_lattice_planes.crystal,
        miller_indices=miller_indices,
        termination_top=termination,
        number_of_repetitions=number_of_layers,
    )

    atomic_layers_repeated_orthogonal_c = AtomicLayersUniqueRepeatedBuilder().get_material(
        atomic_layers_repeated_config
    )

    translation_vector: Vector3D = atomic_layers_repeated_config.get_translation_vector(termination)
    atomic_layers_repeated_terminated = translate(atomic_layers_repeated_orthogonal_c, translation_vector)

    vacuum_configuration = VacuumConfiguration(
        size=vacuum, crystal=atomic_layers_repeated_terminated, direction=AxisEnum.z
    )

    slab_configuration = SlabConfiguration(
        stack_components=[atomic_layers_repeated_terminated, vacuum_configuration],
        direction=AxisEnum.z,
        xy_supercell_matrix=xy_supercell_matrix,
        use_orthogonal_c=use_orthogonal_c,
    )

    builder = SlabBuilder()
    return builder.get_material(slab_configuration)


def create_slab_if_not(material: Material, default_slab_configuration: SlabConfiguration) -> Material:
    slab = material
    if not slab.metadata or slab.metadata["build"]["configuration"]["type"] != SlabConfiguration.__name__:
        print("The material is not a slab. Creating a new slab...")
        slab = create_slab(default_slab_configuration)
    return slab
