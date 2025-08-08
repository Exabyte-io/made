from typing import List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from .types import StackComponentDict
from mat3ra.made.utils import adjust_material_cell_to_set_gap_along_direction
from ....pristine_structures.two_dimensional.slab.helpers import create_slab
from ....pristine_structures.two_dimensional.slab_strained_supercell.builder import SlabStrainedSupercellBuilder
from .....analyze import BaseMaterialAnalyzer
from .....analyze.interface import InterfaceAnalyzer
from .....analyze.slab import SlabMaterialAnalyzer
from .....build_components import MaterialWithBuildMetadata
from .....operations.core.binary import stack


def create_heterostructure(
    stack_component_dicts: List[StackComponentDict],
    gaps: List[float],
    vacuum: float = 10.0,
    use_conventional_cell: bool = True,
    optimize_layer_supercells: bool = True,
) -> MaterialWithBuildMetadata:
    """
    Create a heterostructure by stacking multiple slabs, while applying strain to each slab relative to the first slab.

    Args:
        stack_component_dicts: List of validated stack component configurations
        gaps: List of gaps between adjacent slabs (in Angstroms)
        vacuum: Size of vacuum layer in Angstroms
        use_conventional_cell: Whether to use conventional cell
        optimize_layer_supercells: Whether to find optimal supercells for strained layers

    Returns:
        Heterostructure material with stacked strained slabs
    """
    if len(stack_component_dicts) < 2:
        raise ValueError("At least 2 stack components are required for a heterostructure")

    if len(gaps) != len(stack_component_dicts) - 1:
        raise ValueError("Number of gaps must be one less than number of stack components")

    slabs = []
    for i, component in enumerate(stack_component_dicts):
        slab = create_slab(
            crystal=component.crystal,
            miller_indices=component.miller_indices,
            number_of_layers=component.thickness,
            vacuum=0.0 if i < len(stack_component_dicts) - 1 else vacuum,
            use_conventional_cell=use_conventional_cell,
            xy_supercell_matrix=component.xy_supercell_matrix or [[1, 0], [0, 1]],
        )
        slabs.append(slab)

    strained_slabs = [slabs[0]]  # First slab is the substrate, not strained

    for i in range(1, len(slabs)):
        substrate_slab = strained_slabs[0]
        film_slab = slabs[i]

        substrate_analyzer = SlabMaterialAnalyzer(material=substrate_slab)
        film_analyzer = SlabMaterialAnalyzer(material=film_slab)

        analyzer = InterfaceAnalyzer(
            substrate_slab_configuration=substrate_analyzer.build_configuration,
            film_slab_configuration=film_analyzer.build_configuration,
            substrate_build_parameters=substrate_analyzer.build_parameters,
            film_build_parameters=film_analyzer.build_parameters,
            optimize_film_supercell=optimize_layer_supercells,
        )

        strained_film_config = analyzer.film_strained_configuration

        builder = SlabStrainedSupercellBuilder()
        strained_slab = builder.get_material(strained_film_config)
        strained_slabs.append(strained_slab)

    stacked_materials = []
    for i, slab in enumerate(strained_slabs):
        if i < len(gaps):
            slab_with_gap = adjust_material_cell_to_set_gap_along_direction(slab, gaps[i], AxisEnum.z)
            stacked_materials.append(slab_with_gap)
        else:
            stacked_materials.append(slab)

    heterostructure = stack(stacked_materials, AxisEnum.z)
    heterostructure.name = generate_heterostructure_name(stack_component_dicts)

    return heterostructure


def generate_heterostructure_name(stack_component_dicts: List[StackComponentDict]):
    """Generate a descriptive name for the heterostructure."""

    components = []
    for component in stack_component_dicts:
        analyzer = BaseMaterialAnalyzer(material=component.crystal)
        formula = analyzer.formula
        miller_str = "".join(map(str, component.miller_indices))
        components.append(f"{formula}({miller_str})")

    return f"Heterostructure [{'-'.join(components)}]"
