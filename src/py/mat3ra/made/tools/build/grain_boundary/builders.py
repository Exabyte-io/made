from typing import List

from mat3ra.made.material import Material

from ...analyze import get_chemical_formula
from ..interface import ZSLStrainMatchingInterfaceBuilderParameters, InterfaceConfiguration
from ..interface.builders import ZSLStrainMatchingInterfaceBuilder
from ..supercell import create_supercell
from ..slab import SlabConfiguration, get_terminations, SlabBuilder, SlabSelectorParameters
from .configuration import GrainBoundaryConfiguration


class GrainBoundaryBuilderParameters(ZSLStrainMatchingInterfaceBuilderParameters):
    selected_interface_index: int = 0


class GrainBoundaryBuilder(ZSLStrainMatchingInterfaceBuilder):
    """
    A builder for creating grain boundaries.

    The grain boundary is created by:
    1. creating an interface between two phases,
    2. then rotating the interface by 90 degrees.
    3. Finally, creating a slab from the rotated interface.
    """

    _BuildParametersType: type(GrainBoundaryBuilderParameters) = GrainBoundaryBuilderParameters  # type: ignore
    _ConfigurationType: type(GrainBoundaryConfiguration) = GrainBoundaryConfiguration  # type: ignore
    _GeneratedItemType: type(Material) = Material  # type: ignore
    selector_parameters: type(  # type: ignore
        GrainBoundaryBuilderParameters
    ) = GrainBoundaryBuilderParameters().selected_interface_index  # type: ignore

    def _generate(self, configuration: _ConfigurationType) -> List[Material]:
        interface_config = InterfaceConfiguration(
            film_configuration=configuration.phase_1_configuration,
            substrate_configuration=configuration.phase_2_configuration,
            film_termination=configuration.phase_1_termination,
            substrate_termination=configuration.phase_2_termination,
            distance_z=configuration.gap,
            vacuum=configuration.gap,
        )
        interfaces = super()._generate(interface_config)
        rot_90_degree_matrix = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
        rotated_interfaces = [
            create_supercell(interface, xy_supercell_matrix=rot_90_degree_matrix) for interface in interfaces
        ]
        return rotated_interfaces

    def _finalize(self, materials: List[Material], configuration: _ConfigurationType) -> List[Material]:
        final_slabs = []
        for interface in materials:
            slab_config = SlabConfiguration(
                bulk=interface,
                miller_indices=configuration.slab_configuration.miller_indices,
                thickness=configuration.slab_configuration.thickness,
                vacuum=configuration.slab_configuration.vacuum,
                xy_supercell_matrix=configuration.slab_configuration.xy_supercell_matrix,
                use_conventional_cell=True,
                use_orthogonal_z=True,
            )
            slab_builder = SlabBuilder()
            slab_termination = (
                configuration.slab_termination if configuration.slab_termination else get_terminations(slab_config)[0]
            )
            final_slab = slab_builder.get_material(
                configuration,
                selector_parameters=SlabSelectorParameters(termination=slab_termination),
            )
            final_slabs.append(final_slab)

        materials = super()._finalize(final_slabs, configuration)
        return materials

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        phase_1_formula = get_chemical_formula(configuration.phase_1_configuration.bulk)
        phase_2_formula = get_chemical_formula(configuration.phase_2_configuration.bulk)
        phase_1_miller_indices = "".join([str(i) for i in configuration.phase_1_configuration.miller_indices])
        phase_2_miller_indices = "".join([str(i) for i in configuration.phase_2_configuration.miller_indices])
        new_name = (
            f"{phase_1_formula}({phase_1_miller_indices})-{phase_2_formula}({phase_2_miller_indices}), Grain Boundary"
        )
        material.name = new_name
        return material
