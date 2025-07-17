import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.perturbation import create_perturbation
from mat3ra.made.tools.build.perturbation.builders import PerturbationBuilder
from mat3ra.made.tools.build.perturbation.configuration import PerturbationConfiguration
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.utils.perturbation import SineWavePerturbationFunctionHolder
from mat3ra.utils import assertion as assertion_utils

from .fixtures.monolayer import GRAPHENE


@pytest.mark.parametrize(
    "material_config, supercell_matrix, perturbation_function_params, use_cartesian, coordinates_to_check",
    [
        (
            GRAPHENE,
            [[10, 0, 0], [0, 10, 0], [0, 0, 1]],
            {"amplitude": 0.05, "wavelength": 1},
            False,
            {0: [0.0, 0.0, 0.5], 42: [0.2, 0.1, 0.547552826]},
        ),
    ],
)
def test_sine_perturbation(
    material_config, supercell_matrix, perturbation_function_params, use_cartesian, coordinates_to_check
):
    material = MaterialWithBuildMetadata.create(material_config)
    slab = create_supercell(material, supercell_matrix)

    perturbation_config = PerturbationConfiguration(
        material=slab,
        perturbation_function_holder=SineWavePerturbationFunctionHolder(**perturbation_function_params),
        use_cartesian_coordinates=use_cartesian,
    )
    builder = PerturbationBuilder()
    perturbed_slab = builder.get_material(perturbation_config)
    # Check selected atoms to avoid using 100+ atoms fixture
    for index, expected_coord in coordinates_to_check.items():
        assertion_utils.assert_deep_almost_equal(expected_coord, perturbed_slab.basis.coordinates.values[index])


@pytest.mark.parametrize(
    "material_config, supercell_matrix, perturbation_function_params, use_cartesian, preserve_distance,"
    + " coordinates_to_check, expected_cell",
    [
        (
            GRAPHENE,
            [[10, 0, 0], [0, 10, 0], [0, 0, 1]],
            {"amplitude": 0.05, "wavelength": 1, "phase": 0.25, "axis": "y"},
            False,
            True,
            {0: [0.0, 0.0, 0.5], 42: [0.194051947, 0.1, 0.546942315]},
            [[24.087442, 0.0, 0.0], [-12.043583, 21.367367, 0.0], [0.0, 0.0, 20.0]],
        ),
    ],
)
def test_distance_preserved_sine_perturbation(
    material_config,
    supercell_matrix,
    perturbation_function_params,
    use_cartesian,
    preserve_distance,
    coordinates_to_check,
    expected_cell,
):
    material = MaterialWithBuildMetadata.create(material_config)
    slab = create_supercell(material, supercell_matrix)

    perturbation_config = PerturbationConfiguration(
        material=slab,
        perturbation_function=SineWavePerturbationFunctionHolder(**perturbation_function_params),
        use_cartesian_coordinates=use_cartesian,
    )
    perturbed_slab = create_perturbation(configuration=perturbation_config, preserve_distance=preserve_distance)
    # Check selected atoms to avoid using 100+ atoms fixture
    for index, expected_coord in coordinates_to_check.items():
        assertion_utils.assert_deep_almost_equal(expected_coord, perturbed_slab.basis.coordinates.values[index])
    # Value taken from visually inspected notebook
    assertion_utils.assert_deep_almost_equal(expected_cell, perturbed_slab.lattice.vector_arrays)
