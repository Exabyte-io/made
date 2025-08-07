from typing import List, Tuple, Union

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.utils import adjust_material_cell_to_set_gap_along_direction
from ....pristine_structures.two_dimensional.slab.helpers import create_slab
from ....pristine_structures.two_dimensional.slab_strained_supercell.builder import SlabStrainedSupercellBuilder
from .....analyze import BaseMaterialAnalyzer
from .....analyze.interface import InterfaceAnalyzer
from .....analyze.slab import SlabMaterialAnalyzer
from .....build_components import MaterialWithBuildMetadata
from .....operations.core.binary import stack


def create_heterostructure(
    crystals: List[Union[Material, MaterialWithBuildMetadata]],
    miller_indices: List[Tuple[int, int, int]],
    thicknesses: List[int],
    gaps: List[float],
    vacuum: float = 10.0,
    use_conventional_cell: bool = True,
    optimize_layer_supercells: bool = True,
) -> MaterialWithBuildMetadata:
    """
    Create a heterostructure by stacking multiple slabs, while applying strain to each slab relative to the first slab.
    If `optimize_layer_supercells` is True, it will find larger supercells for strained layers to match the substrate.

    Args:
        crystals: List of crystal materials to create slabs from
        miller_indices: List of Miller indices for each slab surface
        thicknesses: List of layer thicknesses for each slab
        gaps: List of gaps between adjacent slabs (in Angstroms)
        vacuum: Size of vacuum layer in Angstroms
        use_conventional_cell: Whether to use conventional cell
        optimize_layer_supercells: Whether to find optimal supercells for strained layers

    Returns:
        Heterostructure material with stacked strained slabs
    """
    if len(crystals) < 2:
        raise ValueError("At least 2 crystals are required for a heterostructure")

    if len(miller_indices) != len(crystals):
        raise ValueError("Number of Miller indices must match number of crystals")

    if len(thicknesses) != len(crystals):
        raise ValueError("Number of thicknesses must match number of crystals")

    if len(gaps) != len(crystals) - 1:
        raise ValueError("Number of gaps must be one less than number of crystals")

    slabs = []
    for i, (crystal, miller, thickness) in enumerate(zip(crystals, miller_indices, thicknesses)):
        slab = create_slab(
            crystal=crystal,
            miller_indices=miller,
            number_of_layers=thickness,
            vacuum=0.0 if i < len(crystals) - 1 else vacuum,
            use_conventional_cell=use_conventional_cell,
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
    heterostructure.name = generate_heterostructure_name(crystals, miller_indices)

    return heterostructure


def generate_heterostructure_name(crystals, miller_indices):
    """Generate a descriptive name for the heterostructure."""

    components = []
    for crystal, miller in zip(crystals, miller_indices):
        analyzer = BaseMaterialAnalyzer(material=crystal)
        formula = analyzer.formula
        miller_str = "".join(map(str, miller))
        components.append(f"{formula}({miller_str})")

    return f"Heterostructure [{'/'.join(components)}]"
