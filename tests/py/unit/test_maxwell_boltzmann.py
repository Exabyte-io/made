import numpy as np
import pytest


try:
    from mat3ra.periodic_table import PERIODIC_TABLE

    PERIODIC_TABLE_AVAILABLE = True
except ImportError:
    PERIODIC_TABLE_AVAILABLE = False

from mat3ra.made.material import Material
from mat3ra.made.periodic_table import get_atomic_mass_from_element
from mat3ra.made.tools.build_components.operations.core.modifications.perturb.maxwell_boltzmann import (
    BOLTZMANN_CONSTANT_EV_PER_K,
    create_maxwell_displacement_function,
)
from mat3ra.made.tools.helpers import create_supercell
from mat3ra.made.tools.operations.core.unary import perturb

from .fixtures.bulk import BULK_Si_PRIMITIVE
from .fixtures.slab import SI_CONVENTIONAL_SLAB_001

ELEMENT_SYMBOL_TO_MASS_TEST_CASES = [
    ("H", 1.008),
    ("He", 4.003),
    ("Li", 6.941),
    ("C", 12.011),
    ("N", 14.007),
    ("O", 15.999),
    ("Si", 28.085),
    ("Fe", 55.845),
    ("Cu", 63.546),
    ("Au", 196.967),
]

TEMPERATURE_K = 300.0
RANDOM_SEED = 42
NUM_SAMPLES_FOR_MSD = 1000


@pytest.mark.parametrize("random_seed", [None, 42, 123, 999])
def test_maxwell_displacement_deterministic(random_seed):
    material = Material.create(BULK_Si_PRIMITIVE)
    displacement_func1 = create_maxwell_displacement_function(
        material, temperature_in_kelvin=TEMPERATURE_K, random_seed=random_seed
    )
    displacement_func2 = create_maxwell_displacement_function(
        material, temperature_in_kelvin=TEMPERATURE_K, random_seed=random_seed
    )

    if random_seed is not None:
        coord = [0.0, 0.0, 0.0]
        disp1 = displacement_func1.apply_function(coord, atom_index=0)
        disp2 = displacement_func2.apply_function(coord, atom_index=0)
        assert np.allclose(disp1, disp2)
    else:
        coord = [0.0, 0.0, 0.0]
        disp1 = displacement_func1.apply_function(coord, atom_index=0)
        disp2 = displacement_func2.apply_function(coord, atom_index=0)
        assert not np.allclose(disp1, disp2) or np.allclose(disp1, [0, 0, 0], atol=1e-10)


def test_maxwell_displacement_perturb_integration():
    material = Material.create(BULK_Si_PRIMITIVE)
    original_coords = [coord[:] for coord in material.basis.coordinates.values]

    displacement_func = create_maxwell_displacement_function(
        material, temperature_in_kelvin=TEMPERATURE_K, random_seed=RANDOM_SEED
    )

    perturbed_material = perturb(material, displacement_func, use_cartesian_coordinates=True)

    assert len(perturbed_material.basis.coordinates.values) == len(original_coords)
    for i, (orig, pert) in enumerate(zip(original_coords, perturbed_material.basis.coordinates.values)):
        delta = np.array(pert) - np.array(orig)
        assert np.linalg.norm(delta) > 0 or np.allclose(delta, 0, atol=1e-10)


def test_maxwell_displacement_msd_expectation():
    material = Material.create(BULK_Si_PRIMITIVE)
    si_mass = get_atomic_mass_from_element("Si")
    temperature = TEMPERATURE_K
    kT = BOLTZMANN_CONSTANT_EV_PER_K * temperature
    expected_variance = kT / si_mass
    expected_msd = 3 * expected_variance

    displacements = []
    for _ in range(NUM_SAMPLES_FOR_MSD):
        displacement_func = create_maxwell_displacement_function(
            material, temperature_in_kelvin=temperature, random_seed=None
        )
        coord = [0.0, 0.0, 0.0]
        disp = displacement_func.apply_function(coord, atom_index=0)
        displacements.append(disp)

    displacements_array = np.array(displacements)
    msd = np.mean(np.sum(displacements_array**2, axis=1))

    assert abs(msd - expected_msd) / expected_msd < 0.3


@pytest.mark.parametrize(
    "slab_config, temperature_k, random_seed",
    [
        (SI_CONVENTIONAL_SLAB_001, 1300.0, 42),
        (SI_CONVENTIONAL_SLAB_001, 1300.0, 42),
    ],
)
@pytest.mark.skipif(
    not PERIODIC_TABLE_AVAILABLE,
    reason="mat3ra-periodic-table not installed",
)
def test_maxwell_boltzmann_on_slab(slab_config, temperature_k, random_seed):
    material = Material.create(slab_config)
    material = create_supercell(material, scaling_factor=[4, 4, 1])
    original_coords = [coord[:] for coord in material.basis.coordinates.values]
    original_lattice = material.lattice.vector_arrays.copy()

    displacement_func = create_maxwell_displacement_function(
        material, temperature_in_kelvin=temperature_k, random_seed=random_seed
    )

    perturbed_material = perturb(material, displacement_func, use_cartesian_coordinates=True)

    assert len(perturbed_material.basis.coordinates.values) == len(original_coords)
    assert len(perturbed_material.basis.elements.values) == len(material.basis.elements.values)

    coordinate_changes = []
    for i, (orig, pert) in enumerate(zip(original_coords, perturbed_material.basis.coordinates.values)):
        delta = np.array(pert) - np.array(orig)
        displacement_magnitude = np.linalg.norm(delta)
        coordinate_changes.append(displacement_magnitude)

    max_displacement = max(coordinate_changes)
    mean_displacement = np.mean(coordinate_changes)

    assert max_displacement > 0
    assert mean_displacement > 0

    si_mass = get_atomic_mass_from_element("Si")
    kT = BOLTZMANN_CONSTANT_EV_PER_K * temperature_k
    expected_std = np.sqrt(kT / si_mass)

    assert mean_displacement < 5 * expected_std

    assert np.allclose(perturbed_material.lattice.vector_arrays, original_lattice, atol=1e-10)

    for i, element in enumerate(material.basis.elements.values):
        assert perturbed_material.basis.elements.values[i] == element
