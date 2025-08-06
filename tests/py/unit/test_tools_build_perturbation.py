import pytest
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build_components.entities.reusable.three_dimensional.supercell.helpers import create_supercell
from mat3ra.made.tools.build_components.operations.core.modifications.perturb import SineWavePerturbationFunctionHolder
from mat3ra.made.tools.build_components.operations.core.modifications.perturb.builders.base import PerturbationBuilder
from mat3ra.made.tools.build_components.operations.core.modifications.perturb.configuration import (
    PerturbationConfiguration,
)
from mat3ra.made.tools.build_components.operations.core.modifications.perturb.helpers import create_perturbation
from mat3ra.utils import assertion as assertion_utils

from .fixtures.monolayer import GRAPHENE
from .fixtures.nanoribbon.armchair import GRAPHENE_NANORIBBON_ARMCHAIR


@pytest.mark.parametrize(
    "material_config, supercell_matrix, perturbation_function, use_cartesian, coordinates_to_check",
    [
        (
            GRAPHENE,
            [[10, 0, 0], [0, 10, 0], [0, 0, 1]],
            SineWavePerturbationFunctionHolder(amplitude=0.05, wavelength=1),
            False,
            {0: [0.0, 0.0, 0.5], 42: [0.2, 0.1, 0.5475]},
        ),
        (
            GRAPHENE_NANORIBBON_ARMCHAIR,
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            SineWavePerturbationFunctionHolder(amplitude=3.0, wavelength=18.0),
            True,
            {0: [0.1454, 0.279, 0.635], 12: [0.242, 0.3895, 0.6434]},
        ),
    ],
)
def test_perturbation_builder(
    material_config, supercell_matrix, perturbation_function, use_cartesian, coordinates_to_check
):
    material = MaterialWithBuildMetadata.create(material_config)
    slab = create_supercell(material, supercell_matrix)

    perturbation_config = PerturbationConfiguration(
        material=slab,
        perturbation_function_holder=perturbation_function,
        use_cartesian_coordinates=use_cartesian,
    )
    builder = PerturbationBuilder()
    perturbed_slab = builder.get_material(perturbation_config)
    # Check selected atoms to avoid using 100+ atoms fixture
    for index, expected_coord in coordinates_to_check.items():
        assertion_utils.assert_deep_almost_equal(
            expected_coord, perturbed_slab.basis.coordinates.values[index], atol=1e-4
        )


@pytest.mark.parametrize(
    "material_config, supercell_matrix, perturbation_function, use_cartesian, is_isometric,"
    + " coordinates_to_check, expected_cell",
    [
        (
            GRAPHENE,
            [[10, 0, 0], [0, 10, 0], [0, 0, 1]],
            "0.25*sin(x) + 0.25*cos(y)",
            False,
            True,
            {0: [0.0, 0.0, 0.75], 42: [0.19707, 0.10156, 0.7969]},
            [[24.30514, 0.0, 0.0], [-12.15434, 21.03520, 0.0], [0.0, 0.0, 20.0]],
        ),
        (
            GRAPHENE_NANORIBBON_ARMCHAIR,
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "0.3*sin(x*pi)",
            False,
            True,
            {0: [0.1275, 0.2791, 0.5987], 12: [0.2156, 0.3895, 0.6612]},
            [[18.4967, 0.0, 0.0], [0.0, 11.1682, 0.0], [0.0, 0.0, 20.0]],
        ),
    ],
)
def test_create_perturbation(
    material_config,
    supercell_matrix,
    perturbation_function,
    use_cartesian,
    is_isometric,
    coordinates_to_check,
    expected_cell,
):
    material = MaterialWithBuildMetadata.create(material_config)
    slab = create_supercell(material, supercell_matrix)

    perturbed_slab = create_perturbation(
        material=slab,
        perturbation_function=perturbation_function,
        use_cartesian_coordinates=use_cartesian,
        is_isometric=is_isometric,
    )
    # Check selected atoms to avoid using 100+ atoms fixture
    for index, expected_coord in coordinates_to_check.items():
        assertion_utils.assert_deep_almost_equal(
            expected_coord, perturbed_slab.basis.coordinates.values[index], atol=1e-4
        )
    # Value taken from visually inspected notebook
    assertion_utils.assert_deep_almost_equal(expected_cell, perturbed_slab.lattice.vector_arrays, atol=1e-4)
