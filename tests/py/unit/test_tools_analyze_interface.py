from typing import Tuple

import numpy as np
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.utils.assertion import assert_deep_almost_equal

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface import InterfaceAnalyzer
from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.build.slab.builders import AtomicLayersUniqueRepeatedBuilder, SlabBuilder
from mat3ra.made.tools.build.slab.configuration import (
    AtomicLayersUniqueRepeatedConfiguration,
    SlabConfiguration,
    VacuumConfiguration,
)
from mat3ra.made.tools.build.slab.helpers import select_slab_termination, create_slab
from mat3ra.made.tools.operations.core.binary import stack_two_materials
from unit.fixtures.bulk import BULK_Si_CONVENTIONAL, BULK_Si_PRIMITIVE, BULK_Ge_CONVENTIONAL


def get_slab_configuration(
    material_dict: dict,
    miller_indices: Tuple[int, int, int],
    number_of_layers: int,
    vacuum: float,
) -> SlabConfiguration:
    """
    Helper function to create a SlabConfiguration.
    """
    material = Material.model_validate(material_dict)
    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
        material=material, miller_indices=miller_indices
    )
    terminations = crystal_lattice_planes_analyzer.terminations
    termination = select_slab_termination(terminations, None)

    atomic_layers_repeated_configuration = AtomicLayersUniqueRepeatedConfiguration(
        crystal=material,
        miller_indices=miller_indices,
        termination_top=termination,
        number_of_repetitions=number_of_layers,
    )
    atomic_layers_repeated_orthogonal_c = AtomicLayersUniqueRepeatedBuilder().get_material(
        atomic_layers_repeated_configuration
    )
    vacuum_configuration = VacuumConfiguration(
        size=vacuum, crystal=atomic_layers_repeated_orthogonal_c, direction=AxisEnum.z
    )
    return SlabConfiguration(
        stack_components=[atomic_layers_repeated_configuration, vacuum_configuration],
        direction=AxisEnum.z,
    )


def test_interface_analyzer():
    substrate_config = get_slab_configuration(BULK_Si_CONVENTIONAL, (0, 0, 1), 2, 1.0)

    film_config = get_slab_configuration(BULK_Ge_CONVENTIONAL, (0, 0, 1), 2, 5.0)

    analyzer = InterfaceAnalyzer(substrate_config=substrate_config, film_config=film_config)
    substrate_strained_config, film_strained_config = analyzer.get_strained_configurations()

    substrate_slab = SlabBuilder().get_material(substrate_config)
    substrate_strained_slab = SlabBuilder().get_material(substrate_strained_config)
    film_strained_slab = SlabBuilder().get_material(film_strained_config)

    assert_deep_almost_equal(substrate_strained_slab, substrate_slab)

    # Check that lattice vectors orthogonal to the stacking direction are the same
    substrate_in_plane_vectors = substrate_strained_slab.lattice.vector_arrays[:2]
    film_in_plane_vectors = film_strained_slab.lattice.vector_arrays[:2]
    assert_deep_almost_equal(substrate_in_plane_vectors, film_in_plane_vectors)

    stacked_interface = stack_two_materials(
        substrate_strained_slab,
        film_strained_slab,
    )

    assert stacked_interface.basis.elements.values == [
        *substrate_strained_slab.basis.elements.values,
        *film_strained_slab.basis.elements.values,
    ]
