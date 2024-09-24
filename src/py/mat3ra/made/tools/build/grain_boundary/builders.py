from typing import List

import numpy as np
from mat3ra.made.material import Material

from ...analyze import get_chemical_formula
from ..interface import ZSLStrainMatchingInterfaceBuilderParameters, InterfaceConfiguration
from ..interface.builders import ZSLStrainMatchingInterfaceBuilder
from ..supercell import create_supercell
from .configuration import GrainBoundaryConfiguration
from ...modify import add_vacuum
from ...third_party import PymatgenInterface


class GrainBoundaryBuilderParameters(ZSLStrainMatchingInterfaceBuilderParameters):
    default_index: int = 0


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
    _GeneratedItemType: type(PymatgenInterface) = PymatgenInterface  # type: ignore
    selector_parameters: type(  # type: ignore
        GrainBoundaryBuilderParameters
    ) = GrainBoundaryBuilderParameters()  # type: ignore

    def _generate(
        self, configuration: _ConfigurationType
    ) -> List[ZSLStrainMatchingInterfaceBuilder._GeneratedItemType]:
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
        rot_90_degree_matrix = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
        rotated_interfaces = [
            create_supercell(material, supercell_matrix=rot_90_degree_matrix) for material in materials
        ]
        final_slabs: Material = []
        for interface in rotated_interfaces:
            rot_90_degree_matrix = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
            rotated_interfaces = [
                create_supercell(material, supercell_matrix=rot_90_degree_matrix) for material in materials
            ]
            final_slabs = []
            for interface in rotated_interfaces:
                supercell_matrix = np.zeros((3, 3))
                supercell_matrix[:2, :2] = configuration.slab_configuration.xy_supercell_matrix
                supercell_matrix[2, 2] = configuration.slab_configuration.thickness
                final_slab = create_supercell(interface, supercell_matrix=supercell_matrix)
                final_slab_with_vacuum = add_vacuum(final_slab, vacuum=configuration.slab_configuration.vacuum)
                final_slabs.append(final_slab_with_vacuum)

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
