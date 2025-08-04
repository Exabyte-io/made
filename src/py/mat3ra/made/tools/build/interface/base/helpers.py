from typing import List, Optional

from mat3ra.code.array_with_ids import ArrayWithIds

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface import InterfaceAnalyzer
from .builder import InterfaceBuilder
from .build_parameters import InterfaceBuilderParameters
from .configuration import InterfaceConfiguration
from ... import MaterialWithBuildMetadata
from ...vacuum.configuration import VacuumConfiguration
from ....analyze.slab import SlabMaterialAnalyzer


def create_simple_interface_between_slabs(
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
