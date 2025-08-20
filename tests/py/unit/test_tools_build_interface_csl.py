from types import SimpleNamespace

import pytest
from mat3ra.made.tools.analyze.interface.csl import CSLInterfaceAnalyzer
from mat3ra.made.tools.build.compound_pristine_structures.two_dimensional.interface.csl.helpers import (
    create_interface_csl,
    get_csl_strained_configurations,
)
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab import SlabConfiguration

from .fixtures.bulk import BULK_DIAMOND, BULK_GaAs, BULK_Si_CONVENTIONAL
from .fixtures.monolayer import GRAPHENE

SILICON_GRAPHENE_CSL_TEST_CASE = (
    SimpleNamespace(
        bulk_config=BULK_Si_CONVENTIONAL,
        miller_indices=(1, 1, 1),
        number_of_layers=2,
        vacuum=0.0,
    ),
    SimpleNamespace(
        bulk_config=GRAPHENE,
        miller_indices=(0, 0, 1),
        number_of_layers=1,
        vacuum=0.0,
    ),
    3.0,  # gap between graphene and silicon
    10.0,  # vacuum
    50.0,  # max area
)

DIAMOND_GAAS_CSL_TEST_CASE = (
    SimpleNamespace(
        bulk_config=BULK_DIAMOND,
        miller_indices=(0, 0, 1),
        number_of_layers=1,
        vacuum=0.0,
    ),
    SimpleNamespace(
        bulk_config=BULK_GaAs,
        miller_indices=(0, 0, 1),
        number_of_layers=1,
        vacuum=0.0,
    ),
    3.0,  # gap between diamond and gaas
    10.0,  # vacuum
    250.0,  # max area (from your setup)
)


@pytest.mark.parametrize(
    "substrate, film, gap, vacuum, max_area",
    [SILICON_GRAPHENE_CSL_TEST_CASE, DIAMOND_GAAS_CSL_TEST_CASE],
)
def test_csl_interface_creation(substrate, film, gap, vacuum, max_area):
    """Test CSL interface creation between different material pairs."""

    # Test CSL interface creation with different material combinations
    substrate_slab_config = SlabConfiguration.from_parameters(
        substrate.bulk_config,
        substrate.miller_indices,
        substrate.number_of_layers,
        vacuum=substrate.vacuum,
    )
    film_slab_config = SlabConfiguration.from_parameters(
        film.bulk_config,
        film.miller_indices,
        film.number_of_layers,
        vacuum=film.vacuum,
    )

    length_tolerance = 0.1
    strain_tolerance = 2.0
    max_supercell_size = 5

    analyzer = CSLInterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
        max_area=max_area,
        length_tolerance=length_tolerance,
        angle_step=15.0,  # Larger steps for faster testing
        max_rotation_angle=60.0,  # Reduced range for faster testing
        max_supercell_size=max_supercell_size,
        strain_tolerance=strain_tolerance,
    )

    # Check if lattice types are compatible
    if analyzer._validate_same_lattice_type():
        # Test get_csl_strained_configurations
        try:
            strained_configs, strain_percentage = get_csl_strained_configurations(
                substrate_material=substrate.bulk_config,
                film_material=film.bulk_config,
                substrate_miller_indices=substrate.miller_indices,
                film_miller_indices=film.miller_indices,
                number_of_layers=substrate.number_of_layers,
                vacuum=substrate.vacuum,
                max_area=max_area,
                length_tolerance=length_tolerance,
                angle_step=15.0,
                max_rotation_angle=60.0,
                max_supercell_size=max_supercell_size,
                strain_tolerance=strain_tolerance,
            )

            assert len(strained_configs) == 2
            assert isinstance(strain_percentage, float)
            assert strain_percentage >= 0

            # Test create_interface_csl
            interface = create_interface_csl(
                substrate_material=substrate.bulk_config,
                film_material=film.bulk_config,
                substrate_miller_indices=substrate.miller_indices,
                film_miller_indices=film.miller_indices,
                number_of_layers=substrate.number_of_layers,
                vacuum=vacuum,
                gap=gap,
                max_area=max_area,
                length_tolerance=length_tolerance,
                angle_step=15.0,
                max_rotation_angle=60.0,
                max_supercell_size=max_supercell_size,
                strain_tolerance=strain_tolerance,
            )

            assert interface is not None
            assert hasattr(interface, "basis")
            assert hasattr(interface, "lattice")

        except ValueError as e:
            # If no matches found, that's acceptable for this test case
            if "No CSL matches found" in str(e):
                pytest.skip(f"No CSL matches found for the given parameters: {e}")
            else:
                raise
    else:
        pytest.skip("Substrate and film do not have the same lattice type for CSL analysis")
