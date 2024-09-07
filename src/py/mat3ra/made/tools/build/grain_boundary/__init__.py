from pydantic import BaseModel
from mat3ra.made.material import Material
from mat3ra.code.entity import InMemoryEntity

from .. import BaseBuilder
from ..interface import ZSLStrainMatchingInterfaceBuilderParameters, InterfaceConfiguration
from ..supercell import create_supercell
from ...analyze import get_chemical_formula
from ..interface.builders import ZSLStrainMatchingInterfaceBuilder
from ..slab import SlabConfiguration, Termination, get_terminations, SlabBuilder, SlabSelectorParameters


class GrainBoundaryConfiguration(BaseModel, InMemoryEntity):
    phase_1_configuration: SlabConfiguration
    phase_2_configuration: SlabConfiguration
    phase_1_termination: Termination
    phase_2_termination: Termination
    gap: float = 3.0
    slab_configuration: SlabConfiguration
    slab_termination: Termination

    @property
    def _json(self):
        return {
            "type": "GrainBoundaryConfiguration",
            "phase_1_configuration": self.phase_1_configuration.to_json(),
            "phase_2_configuration": self.phase_2_configuration.to_json(),
            "phase_1_termination": str(self.phase_1_termination),
            "phase_2_termination": str(self.phase_2_termination),
            "gap": self.gap,
            "slab_configuration": self.slab_configuration.to_json(),
        }


class GrainBoundaryBuilderParameters(ZSLStrainMatchingInterfaceBuilderParameters):
    selected_interface_index: int = 0


class GrainBoundaryBuilder(BaseBuilder):
    _BuildParametersType = GrainBoundaryBuilderParameters
    _ConfigurationType: type(GrainBoundaryConfiguration) = GrainBoundaryConfiguration  # type: ignore

    def _generate(self, configuration: GrainBoundaryConfiguration):
        interface = self._create_interface(configuration)
        rotated_interface = create_supercell(interface, [[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
        slab_config = SlabConfiguration(
            bulk=rotated_interface,
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

        return [final_slab]

    def _create_interface(self, configuration: GrainBoundaryConfiguration):
        phase_1_config = configuration.phase_1_configuration
        phase_2_config = configuration.phase_2_configuration
        phase_1_termination = (
            configuration.phase_1_termination
            if configuration.phase_1_termination
            else get_terminations(phase_1_config)[0]
        )
        phase_2_termination = (
            configuration.phase_2_termination
            if configuration.phase_2_termination
            else get_terminations(phase_2_config)[0]
        )

        interface_configuration = InterfaceConfiguration(
            film_configuration=phase_1_config,
            substrate_configuration=phase_2_config,
            film_termination=phase_1_termination,
            substrate_termination=phase_2_termination,
            distance_z=configuration.gap,
            vacuum=configuration.gap,
        )

        zsl_parameters = self.build_parameters
        matched_interfaces_builder = ZSLStrainMatchingInterfaceBuilder(
            build_parameters=ZSLStrainMatchingInterfaceBuilderParameters(strain_matching_parameters=zsl_parameters)
        )

        interfaces_sorted_by_size_and_strain = matched_interfaces_builder.get_materials(
            configuration=interface_configuration
        )

        selected_interface = interfaces_sorted_by_size_and_strain[self.build_parameters.selected_interface_index]
        return selected_interface

    def _update_material_name(self, material: Material, configuration: GrainBoundaryConfiguration) -> Material:
        phase_1_formula = get_chemical_formula(configuration.phase_1_configuration.bulk)
        phase_2_formula = get_chemical_formula(configuration.phase_2_configuration.bulk)
        phase_1_miller_indices = "".join([str(i) for i in configuration.phase_1_configuration.miller_indices])
        phase_2_miller_indices = "".join([str(i) for i in configuration.phase_2_configuration.miller_indices])
        new_name = (
            f"{phase_1_formula}({phase_1_miller_indices})-{phase_2_formula}({phase_2_miller_indices}), Grain Boundary"
        )
        material.name = new_name
        return material
