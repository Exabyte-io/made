from typing import List, Optional, Union

from mat3ra.made.material import Material
from .builders import SlabBuilder, SlabSelectorParameters, SlabBuilderParameters
from .configuration import (
    SlabConfiguration, 
    AtomicLayersUnique, 
    VacuumConfiguration,
    MillerIndicesSchema
)
from .termination import Termination

CACHED_BUILDER = None


def get_terminations(
    configuration: SlabConfiguration, build_parameters: Optional[SlabBuilderParameters] = None
) -> List[Termination]:
    global CACHED_BUILDER
    CACHED_BUILDER = SlabBuilder(build_parameters=build_parameters)
    return CACHED_BUILDER.get_terminations(configuration)


def create_slab(
    configuration: Optional[SlabConfiguration] = None,
    build_parameters: Optional[SlabBuilderParameters] = None,
    use_cached_builder: bool = False,
    # Individual parameters for convenience
    crystal: Optional[Material] = None,
    miller_indices: Optional[Union[MillerIndicesSchema, List[int], tuple]] = None,
    use_conventional_cell: Optional[bool] = None,
    termination: Optional[Termination] = None,
    number_of_layers: Optional[int] = None,
    vacuum: Optional[float] = None,
    xy_supercell_matrix: Optional[List[List[int]]] = None,
) -> Material:
    # If configuration is provided, use it directly
    if configuration is not None:
        builder = (
            CACHED_BUILDER if use_cached_builder and CACHED_BUILDER else SlabBuilder(build_parameters=build_parameters)
        )
        return builder.get_material(configuration)
    
    # Otherwise, construct configuration from individual parameters
    if crystal is None or miller_indices is None:
        raise ValueError("Either 'configuration' or both 'crystal' and 'miller_indices' must be provided")
    
    # Convert miller_indices to MillerIndicesSchema if needed
    if isinstance(miller_indices, (list, tuple)):
        miller_indices_schema = MillerIndicesSchema(root=list(miller_indices))
    else:
        miller_indices_schema = miller_indices
    
    # Create atomic layers configuration
    atomic_layers = AtomicLayersUnique(
        crystal=crystal,
        miller_indices=miller_indices_schema,
        use_conventional_cell=use_conventional_cell or False,
    )
    
    # Create vacuum configuration
    vacuum_config = VacuumConfiguration(
        size=vacuum or 0.0,
        crystal=crystal
    )
    
    # Create slab configuration
    slab_config = SlabConfiguration(
        stack_components=[atomic_layers, vacuum_config],
        xy_supercell_matrix=xy_supercell_matrix or [[1, 0], [0, 1]],
    )
    
    builder = SlabBuilder(build_parameters=build_parameters)
    return builder.get_material(slab_config)


def create_slab_if_not(material: Material, default_slab_configuration: SlabConfiguration) -> Material:
    slab = material
    if not slab.metadata or slab.metadata["build"]["configuration"]["type"] != SlabConfiguration.__name__:
        print("The material is not a slab. Creating a new slab...")
        slab = create_slab(default_slab_configuration)
    return slab
