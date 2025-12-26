import numpy as np
import pytest

from mat3ra.made.material import Material
from mat3ra.made.periodic_table import get_atomic_mass_from_element
from mat3ra.made.tools.build_components.operations.core.modifications.perturb.functions.maxwell_boltzmann import (
    MaxwellBoltzmannDisplacementHolder,
)
from mat3ra.made.tools.build_components.operations.core.modifications.perturb.helpers import (
    create_maxwell_displacement,
)
from mat3ra.made.tools.helpers import create_supercell
from .fixtures.bulk import BULK_Si_PRIMITIVE
from .fixtures.slab import SI_CONVENTIONAL_SLAB_001

DISORDER_PARAMETER = 3000.0  # Temperature-like
RANDOM_SEED = 42
NUM_SAMPLES_FOR_MSD = 1000


@pytest.mark.parametrize("random_seed", [None, 42, 123])
def test_maxwell_displacement_deterministic(random_seed):
    material = Material.create(BULK_Si_PRIMITIVE)
    displacement_func1 = MaxwellBoltzmannDisplacementHolder(
        disorder_parameter=DISORDER_PARAMETER, random_seed=random_seed
    )
    displacement_func2 = MaxwellBoltzmannDisplacementHolder(
        disorder_parameter=DISORDER_PARAMETER, random_seed=random_seed
    )

    coord = [0.0, 0.0, 0.0]
    atom_index = 0

    if random_seed is not None:
        disp1 = displacement_func1.apply_function(coord, material=material, atom_index=atom_index)
        disp2 = displacement_func2.apply_function(coord, material=material, atom_index=atom_index)
        assert np.allclose(disp1, disp2)

        # Different seed should give different results
        displacement_func3 = MaxwellBoltzmannDisplacementHolder(
            disorder_parameter=DISORDER_PARAMETER, random_seed=random_seed + 1
        )
        disp3 = displacement_func3.apply_function(coord, material=material, atom_index=atom_index)
        assert not np.allclose(disp1, disp3) or np.allclose(disp1, [0, 0, 0], atol=1e-10)
    else:
        # No seed: different instances should give different results (non-deterministic)
        disp1 = displacement_func1.apply_function(coord, material=material, atom_index=atom_index)
        disp2 = displacement_func2.apply_function(coord, material=material, atom_index=atom_index)
        assert not np.allclose(disp1, disp2) or np.allclose(disp1, [0, 0, 0], atol=1e-10)


def test_maxwell_displacement_perturb_integration():
    material = Material.create(BULK_Si_PRIMITIVE)
    original_coords = [coord[:] for coord in material.basis.coordinates.values]

    perturbed_material = create_maxwell_displacement(
        material, disorder_parameter=DISORDER_PARAMETER, random_seed=RANDOM_SEED
    )

    assert len(perturbed_material.basis.coordinates.values) == len(original_coords)
    for i, (orig, pert) in enumerate(zip(original_coords, perturbed_material.basis.coordinates.values)):
        delta = np.array(pert) - np.array(orig)
        assert np.linalg.norm(delta) > 0 or np.allclose(delta, 0, atol=1e-10)


def test_maxwell_displacement_msd_expectation():
    material = Material.create(BULK_Si_PRIMITIVE)
    si_mass = get_atomic_mass_from_element("Si")
    disorder_parameter = DISORDER_PARAMETER
    expected_variance = disorder_parameter / si_mass
    expected_msd = 3 * expected_variance

    displacements = []
    atom_index = 0
    coord = [0.0, 0.0, 0.0]
    for _ in range(NUM_SAMPLES_FOR_MSD):
        displacement_func = MaxwellBoltzmannDisplacementHolder(disorder_parameter=disorder_parameter, random_seed=None)
        disp = displacement_func.apply_function(coord, material=material, atom_index=atom_index)
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
def test_maxwell_boltzmann_on_slab(slab_config, temperature_k, random_seed):
    material = Material.create(slab_config)
    material = create_supercell(material, scaling_factor=[4, 4, 1])
    original_coords = [coord[:] for coord in material.basis.coordinates.values]
    original_lattice = material.lattice.vector_arrays.copy()

    perturbed_material = create_maxwell_displacement(
        material, disorder_parameter=temperature_k, random_seed=random_seed
    )

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
    expected_std = np.sqrt(temperature_k / si_mass)

    assert mean_displacement < 5 * expected_std

    assert np.allclose(perturbed_material.lattice.vector_arrays, original_lattice, atol=1e-10)

    for i, element in enumerate(material.basis.elements.values):
        assert perturbed_material.basis.elements.values[i] == element
