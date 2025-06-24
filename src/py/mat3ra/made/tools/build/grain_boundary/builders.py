from typing import List

import numpy as np

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.commensurate import CommensurateLatticeInterfaceAnalyzer
from .configuration import SurfaceGrainBoundaryConfiguration, SlabGrainBoundaryConfiguration
from ..interface import ZSLStrainMatchingInterfaceBuilderParameters, InterfaceConfiguration
from ..interface.builders import (
    ZSLStrainMatchingInterfaceBuilder,
    InterfaceBuilder,
    CommensurateLatticeInterfaceBuilderParameters,
)
from ..slab.configurations import SlabConfiguration
from ..slab.builders import SlabBuilder, SlabBuilderParameters
from ..supercell import create_supercell
from ..utils import stack_two_materials_xy
from ...analyze.other import get_chemical_formula
from ...third_party import PymatgenInterface


class SlabGrainBoundaryBuilderParameters(ZSLStrainMatchingInterfaceBuilderParameters):
    default_index: int = 0


class SlabGrainBoundaryBuilder(ZSLStrainMatchingInterfaceBuilder):
    """
    A builder for creating grain boundaries.

    The grain boundary is created by:
    1. creating an interface between two phases,
    2. then rotating the interface by 90 degrees.
    3. Finally, creating a slab from the rotated interface.
    """

    _BuildParametersType: type(SlabGrainBoundaryBuilderParameters) = SlabGrainBoundaryBuilderParameters  # type: ignore
    _DefaultBuildParameters = SlabGrainBoundaryBuilderParameters()
    _ConfigurationType: type(SlabGrainBoundaryConfiguration) = SlabGrainBoundaryConfiguration  # type: ignore
    _GeneratedItemType: PymatgenInterface = PymatgenInterface  # type: ignore
    selector_parameters: type(  # type: ignore
        SlabGrainBoundaryBuilderParameters
    ) = SlabGrainBoundaryBuilderParameters()  # type: ignore

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:  # type: ignore
        interface_config = InterfaceConfiguration(
            film_configuration=configuration.phase_1_configuration,
            substrate_configuration=configuration.phase_2_configuration,
            film_termination=configuration.phase_1_termination,
            substrate_termination=configuration.phase_2_termination,
            distance_z=configuration.gap,
            vacuum=configuration.gap,
        )
        return super()._generate(interface_config)

    def _finalize(self, materials: List[Material], configuration: _ConfigurationType) -> List[Material]:
        rot_90_degree_matrix = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
        rotated_interfaces = [
            create_supercell(material, supercell_matrix=rot_90_degree_matrix) for material in materials
        ]
        final_slabs: List[Material] = []
        for interface in rotated_interfaces:
            final_slab_config = SlabConfiguration.from_parameters(
                material_or_dict=interface,
                miller_indices=configuration.slab_configuration.miller_indices,
                number_of_layers=configuration.slab_configuration.number_of_layers,
                vacuum=configuration.slab_configuration.vacuum,
            )
            slab_builder_parameters = SlabBuilderParameters(
                xy_supercell_matrix=configuration.slab_configuration.xy_supercell_matrix,
                use_orthogonal_c=True,
            )
            builder = SlabBuilder(build_parameters=slab_builder_parameters)
            final_slab = builder.get_material(final_slab_config)
            final_slabs.append(final_slab)

        return super()._finalize(final_slabs, configuration)

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


class SurfaceGrainBoundaryBuilderParameters(CommensurateLatticeInterfaceBuilderParameters):
    """
    Parameters for creating a grain boundary between two surface phases.

    Args:
        edge_inclusion_tolerance (float): The tolerance to include atoms on the edge of each phase, in angstroms.
        distance_tolerance (float): The distance tolerance to remove atoms that are too close, in angstroms.
    """

    edge_inclusion_tolerance: float = 1.0
    distance_tolerance: float = 1.0


class SurfaceGrainBoundaryBuilder(InterfaceBuilder):
    _ConfigurationType: type(SurfaceGrainBoundaryConfiguration) = SurfaceGrainBoundaryConfiguration  # type: ignore
    _BuildParametersType = SurfaceGrainBoundaryBuilderParameters
    _DefaultBuildParameters = SurfaceGrainBoundaryBuilderParameters()

    def _generate(self, configuration: SurfaceGrainBoundaryConfiguration) -> Material:
        # Create slab configuration from the material
        slab_config = SlabConfiguration.from_parameters(
            configuration.film, miller_indices=(0, 0, 1), number_of_layers=1, vacuum=0.0
        )

        analyzer = CommensurateLatticeInterfaceAnalyzer(
            substrate_slab_configuration=slab_config,
            target_angle=configuration.twist_angle,
            angle_tolerance=self.build_parameters.angle_tolerance,
        )

        match_holder = analyzer.commensurate_lattice_match_holders[0]
        matrix1 = np.array(match_holder.xy_supercell_matrix_substrate)
        matrix2 = np.array(match_holder.xy_supercell_matrix_film)

        if configuration.xy_supercell_matrix != [[1, 0], [0, 1]]:
            base_matrix = np.array(configuration.xy_supercell_matrix)
            matrix1 = np.dot(base_matrix, matrix1)
            matrix2 = np.dot(base_matrix, matrix2)

        phase_1_material = create_supercell(configuration.film, matrix1.tolist())
        phase_2_material = create_supercell(configuration.film, matrix2.tolist())

        grain_boundary = stack_two_materials_xy(
            phase_1_material,
            phase_2_material,
            gap=configuration.gap,
            edge_inclusion_tolerance=self.build_parameters.edge_inclusion_tolerance,
            distance_tolerance=self.build_parameters.distance_tolerance,
        )

        return grain_boundary

    def _update_material_name(self, material: Material, configuration: SurfaceGrainBoundaryConfiguration) -> Material:
        new_name = f"{configuration.film.name}, Grain Boundary ({configuration.twist_angle:.2f}°)"
        material.name = new_name
        return material
