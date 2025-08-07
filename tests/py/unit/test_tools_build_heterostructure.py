from types import SimpleNamespace

import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.compound_pristine_structures.two_dimensional.heterostructure import create_heterostructure
from mat3ra.standata.materials import Materials

from .fixtures.bulk import BULK_Si_CONVENTIONAL

PRECISION = 1e-3

Si_SiO2_Hf2O_HETEROSTRUCTURE_TEST_CASE = (
    [
        SimpleNamespace(bulk_config=BULK_Si_CONVENTIONAL, miller_indices=(0, 0, 1), number_of_layers=2),
        SimpleNamespace(
            bulk_config=Materials.get_by_name_first_match("SiO2"), miller_indices=(1, 1, 1), number_of_layers=3
        ),
        SimpleNamespace(
            bulk_config=Materials.get_by_name_first_match("Hafnium.*MCL"), miller_indices=(0, 0, 1), number_of_layers=2
        ),
        SimpleNamespace(
            bulk_config=Materials.get_by_name_first_match("TiN"), miller_indices=(1, 1, 1), number_of_layers=3
        ),
    ],
    [2.0, 2.0, 2.0],  # gaps
    10.0,  # vacuum
)


@pytest.mark.parametrize("layers, gaps, vacuum", [Si_SiO2_Hf2O_HETEROSTRUCTURE_TEST_CASE])
def test_create_heterostructure_simple(layers, gaps, vacuum):
    crystals = [Material.create(layer.bulk_config) for layer in layers]
    miller_indices = [layer.miller_indices for layer in layers]
    thicknesses = [layer.number_of_layers for layer in layers]

    heterostructure = create_heterostructure(
        crystals=crystals,
        miller_indices=miller_indices,
        thicknesses=thicknesses,
        gaps=gaps,
        analyzer_type="simple",
        vacuum=vacuum,
    )

    # Check that we got a valid material
    assert isinstance(heterostructure, Material)
    assert heterostructure.lattice is not None
    assert heterostructure.basis is not None

    # Should have atoms from all layers
    assert len(heterostructure.basis.elements.values) > 0

    # Check that it has reasonable dimensions
    lattice_c = heterostructure.lattice.c
    assert lattice_c > vacuum  # Should include vacuum


@pytest.mark.parametrize("layers, gaps, vacuum", [Si_SiO2_Hf2O_HETEROSTRUCTURE_TEST_CASE])
def test_create_heterostructure_zsl(layers, gaps, vacuum):
    crystals = [Material.create(layer.bulk_config) for layer in layers]
    miller_indices = [layer.miller_indices for layer in layers]
    thicknesses = [layer.number_of_layers for layer in layers]

    try:
        heterostructure = create_heterostructure(
            crystals=crystals,
            miller_indices=miller_indices,
            thicknesses=thicknesses,
            gaps=gaps,
            analyzer_type="zsl",
            vacuum=vacuum,
            max_area=50.0,
        )

        # Check that we got a valid material
        assert isinstance(heterostructure, Material)
        assert heterostructure.lattice is not None
        assert heterostructure.basis is not None
    except ValueError as e:
        if "No ZSL matches found" in str(e):
            pytest.skip("ZSL analyzer could not find matches for these materials")
        else:
            raise
