from mat3ra.made.material import Material
from sympy import ceiling

from ...slab.helpers import create_slab

from ....analyze.slab import SlabMaterialAnalyzer
from ...defect.configuration import SlabDefectConfigurationLegacy
from ...defect.slab.builders import SlabDefectBuilder
from ...defect.slab.configuration import SlabDefectConfiguration
from ...slab.builders import SlabBuilder

from src.py.mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from ....modify import filter_by_box


def create_slab_defect(slab: Material, isolated_defect: Material, vacuum_thickness: float = 0.0) -> Material:
    """
    Create a slab defect by merging an isolated defect with a slab material.
    Args:
        slab: The original slab material.
        isolated_defect: The isolated defect material.
        additional_layers: Number of additional layers to add to the slab for defect creation.
    Returns:
        Material: The new slab material with additional layers and vacuum if needed.
    """

    analyzer = SlabMaterialAnalyzer(material=slab)

    vacuum_configuration = VacuumConfiguration(vacuum_thickness=vacuum_thickness)

    return stack(slab, isolated_defect, vacuum_configuration)


def create_slab_with_additional_layers(
    slab: Material, number_of_additional_layers: int = 1, vacuum: float = 0.0
) -> Material:
    """
    Create a slab with additional layers by merging a slab material with additional layers and vacuum.
    Args:
        slab: The original slab material.
        number_of_additional_layers: Number of additional layers to add to the slab.
        vacuum: Thickness of the vacuum layer to be added.
    Returns:
        Material: The new slab material with additional layers and vacuum if needed.
    """
    # get_slab_build_configuration handles the supercell creation, etc.
    slab_configuration = get_slab_build_configuration(slab.metadata)
    total_number_of_layers = slab_configuration.number_of_layers + number_of_additional_layers
    return create_slab(
        crystal=slab,
        miller_indices=slab_configuration.miller_indices,
        termination=slab_configuration.termination,
        number_of_layers=total_number_of_layers,
        vacuum=vacuum,
        xy_supercell_matrix=slab_configuration.xy_supercell_matrix,
    )


def create_slab_with_additional_layers_float(
    slab: Material, number_of_additional_layers: float = 0.5, vacuum: float = 0.0
) -> Material:
    """
    Create a slab with additional layers by merging a slab material with additional layers and vacuum.
    Args:
        slab: The original slab material.
        number_of_additional_layers: Number of additional layers to add to the slab.
        vacuum: Thickness of the vacuum layer to be added.
    Returns:
        Material: The new slab material with additional layers and vacuum if needed.
    """
    slab_analyzer = SlabMaterialAnalyzer(material=slab)
    ceiling_number_of_additional_layers = int(ceiling(number_of_additional_layers))
    slab_with_int_additional_layers = create_slab_with_additional_layers(
        slab=slab,
        number_of_additional_layers=ceiling_number_of_additional_layers,
        vacuum=vacuum,
    )
    vacuum_thickness_in_layers = slab_analyzer.vacuum_thickness_in_layers
    number_of_layers = slab_analyzer.number_of_layers
    max_z_crystal_coordinate = (number_of_additional_layers + number_of_layers + vacuum_thickness_in_layers) / (
        ceiling_number_of_additional_layers + number_of_layers + vacuum_thickness_in_layers
    )
    return filter_by_box(
        slab_with_int_additional_layers,
        max_coordinate=[1, 1, max_z_crystal_coordinate],
        reset_ids=True,
    )
