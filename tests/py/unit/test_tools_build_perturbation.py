from mat3ra.made.material import Material
from mat3ra.made.tools.build.perturbation import create_perturbation
from mat3ra.made.tools.build.perturbation.builders import SlabPerturbationBuilder
from mat3ra.made.tools.build.perturbation.configuration import PerturbationConfiguration
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.utils.perturbation import SineWavePerturbationFunctionHolder
from mat3ra.utils import assertion as assertion_utils

from .fixtures.monolayer import GRAPHENE


def test_sine_perturbation():
    material = Material.create(GRAPHENE)
    slab = create_supercell(material, [[10, 0, 0], [0, 10, 0], [0, 0, 1]])

    perturbation_config = PerturbationConfiguration(
        material=slab,
        perturbation_function=SineWavePerturbationFunctionHolder(amplitude=0.05, wavelength=1),
        use_cartesian_coordinates=False,
    )
    builder = SlabPerturbationBuilder()
    perturbed_slab = builder.get_material(perturbation_config)
    # Check selected atoms to avoid using 100+ atoms fixture
    assertion_utils.assert_deep_almost_equal([0.0, 0.0, 0.5], perturbed_slab.basis.coordinates.values[0])
    assertion_utils.assert_deep_almost_equal([0.2, 0.1, 0.547552826], perturbed_slab.basis.coordinates.values[42])


def test_distance_preserved_sine_perturbation():
    material = Material.create(GRAPHENE)
    slab = create_supercell(material, [[10, 0, 0], [0, 10, 0], [0, 0, 1]])

    perturbation_config = PerturbationConfiguration(
        material=slab,
        perturbation_function=SineWavePerturbationFunctionHolder(amplitude=0.05, wavelength=1, phase=0.25, axis="y"),
        use_cartesian_coordinates=False,
    )
    perturbed_slab = create_perturbation(configuration=perturbation_config, preserve_distance=True)
    # Check selected atoms to avoid using 100+ atoms fixture
    assertion_utils.assert_deep_almost_equal([0.0, 0.0, 0.5], perturbed_slab.basis.coordinates.values[0])
    assertion_utils.assert_deep_almost_equal(
        [0.194051947, 0.1, 0.546942315], perturbed_slab.basis.coordinates.values[42]
    )
    # Value taken from visually inspected notebook
    expected_cell = [
        [24.087442, 0.0, 0.0],
        [-12.043583, 21.367367, 0.0],
        [0.0, 0.0, 20.0],
    ]
    assertion_utils.assert_deep_almost_equal(expected_cell, perturbed_slab.lattice.vector_arrays)
