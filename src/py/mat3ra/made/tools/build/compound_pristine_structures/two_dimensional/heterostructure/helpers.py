from typing import List

import numpy as np
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from .types import StackComponentDict
from mat3ra.made.utils import adjust_material_cell_to_set_gap_along_direction
from ....pristine_structures.two_dimensional.slab.helpers import create_slab
from ....pristine_structures.two_dimensional.slab_strained_supercell.builder import SlabStrainedSupercellBuilder
from ....pristine_structures.two_dimensional.slab_strained_supercell.configuration import SlabStrainedSupercellConfiguration
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

    # Process all materials to find a common supercell that works for all
    # We'll use the first material as the reference substrate
    reference_substrate = slabs[0]
    strained_slabs = []
    common_supercell_matrix = None

    for i in range(1, len(slabs)):
        film_slab = slabs[i]

        substrate_analyzer = SlabMaterialAnalyzer(material=reference_substrate)
        film_analyzer = SlabMaterialAnalyzer(material=film_slab)

        analyzer = InterfaceAnalyzer(
            substrate_slab_configuration=substrate_analyzer.build_configuration,
            film_slab_configuration=film_analyzer.build_configuration,
            substrate_build_parameters=substrate_analyzer.build_parameters,
            film_build_parameters=film_analyzer.build_parameters,
            optimize_film_supercell=optimize_layer_supercells,
        )

        film_config = analyzer.film_strained_configuration
        if common_supercell_matrix is None:
            common_supercell_matrix = film_config.xy_supercell_matrix

    # Second pass: create all materials with the common supercell
    substrate_analyzer = SlabMaterialAnalyzer(material=reference_substrate)

    # Create substrate with common supercell
    substrate_config = SlabStrainedSupercellConfiguration(
        **substrate_analyzer.build_configuration.model_dump(),
        strain_matrix=Matrix3x3Schema(root=np.eye(3).tolist()),
        xy_supercell_matrix=common_supercell_matrix,
    )

    builder = SlabStrainedSupercellBuilder()
    strained_substrate = builder.get_material(substrate_config)
    strained_slabs.append(strained_substrate)

    # Create all films with strain to match the substrate
    for i in range(1, len(slabs)):
        film_slab = slabs[i]

        substrate_analyzer = SlabMaterialAnalyzer(material=reference_substrate)
        film_analyzer = SlabMaterialAnalyzer(material=film_slab)

        analyzer = InterfaceAnalyzer(
            substrate_slab_configuration=substrate_analyzer.build_configuration,
            film_slab_configuration=film_analyzer.build_configuration,
            substrate_build_parameters=substrate_analyzer.build_parameters,
            film_build_parameters=film_analyzer.build_parameters,
            optimize_film_supercell=optimize_layer_supercells,
        )

        film_config = analyzer.film_strained_configuration

        # Override the supercell matrix to use the common one
        film_config_dict = film_config.model_dump()
        film_config_dict['xy_supercell_matrix'] = common_supercell_matrix
        film_config_common = SlabStrainedSupercellConfiguration(**film_config_dict)

        strained_film = builder.get_material(film_config_common)
        strained_slabs.append(strained_film)

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
