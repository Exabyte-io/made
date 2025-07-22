import pytest
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.adatom.helpers import create_adatom_defect, create_multiple_adatom_defects
from mat3ra.made.tools.build.defect.enums import AdatomPlacementMethodEnum
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001


@pytest.mark.parametrize(
    "crystal_config, position_on_surface, distance_z, chemical_element, adatom_placement_method, expected_last_coord",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            [0.5, 0.5],
            2.0,
            "Si",
            AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE,
            [0.5, 0.5, 0.5748],
        ),
        (
            SI_CONVENTIONAL_SLAB_001,
            [0.55, 0.51],
            2.0,
            "Si",
            AdatomPlacementMethodEnum.EXACT_COORDINATE,
            [0.55, 0.51, 0.6231],
        ),
    ],
)
def test_create_adatom(
    crystal_config,
    position_on_surface,
    distance_z,
    chemical_element,
    adatom_placement_method,
    expected_last_coord,
):
    slab = MaterialWithBuildMetadata.create(crystal_config)
    defect = create_adatom_defect(slab, position_on_surface, distance_z, adatom_placement_method, chemical_element)

    assertion_utils.assert_deep_almost_equal(expected_last_coord, defect.basis.coordinates.values[-1], atol=1e-4)


@pytest.mark.parametrize(
    "crystal_config, adatom_dicts, placement_method, expected_num_atoms",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            [
                {"element": "Si", "coordinate": [0.5, 0.5], "distance_z": 2.0},
                {"element": "C", "coordinate": [0.25, 0.25], "distance_z": 1.5},
                {"element": "N", "coordinate": [0.75, 0.75]},  # distance_z will default to 1.0
            ],
            AdatomPlacementMethodEnum.EXACT_COORDINATE,
            19,  # 16 atoms in slab + 3 adatoms
        ),
    ],
)
def test_create_multiple_adatom_defects(crystal_config, adatom_dicts, placement_method, expected_num_atoms):
    slab = MaterialWithBuildMetadata.create(crystal_config)
    defects = create_multiple_adatom_defects(slab, adatom_dicts, placement_method)

    assert len(defects.basis.elements.values) == expected_num_atoms
    
    # Check that all expected elements are present
    expected_elements = ["Si"] * 16 + ["Si", "C", "N"]  # Original slab + added adatoms
    actual_elements = defects.basis.elements.values
    assert sorted(actual_elements) == sorted(expected_elements)
