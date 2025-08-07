import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defective_structures.two_dimensional.adatom.helpers import (
    create_defect_adatom,
    create_multiple_adatom_defects,
)
from mat3ra.made.tools.build.defective_structures.two_dimensional.adatom.types import AdatomDefectDict
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.helpers import create_slab
from mat3ra.made.tools.build_components.operations.core.combinations.enums import AdatomPlacementMethodEnum
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.adatom import SLAB_Si_3_ADATOMS
from unit.fixtures.bulk import BULK_Si_CONVENTIONAL
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001
from unit.utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "crystal_config, position_on_surface, distance_z, chemical_element, adatom_placement_method, expected_last_coord",
    [
        (
            BULK_Si_CONVENTIONAL,
            [0.5, 0.5],
            2.0,
            "Si",
            AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE.value,
            [0.5, 0.5, 0.4659],
        ),
        (
            BULK_Si_CONVENTIONAL,
            [0.55, 0.51],
            2.0,
            "Si",
            AdatomPlacementMethodEnum.EXACT_COORDINATE.value,
            [0.55, 0.51, 0.4621],
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
    slab = create_slab(
        crystal=Material.create(crystal_config),
        number_of_layers=2,
        xy_supercell_matrix=[[2, 0], [0, 2]],
        termination_top_formula=None,
        termination_bottom_formula=None,
    )
    defect = create_defect_adatom(slab, position_on_surface, distance_z, adatom_placement_method, chemical_element)

    assertion_utils.assert_deep_almost_equal(expected_last_coord, defect.basis.coordinates.values[-1], atol=1e-4)


@pytest.mark.parametrize(
    "crystal_config, adatom_dicts, placement_method, expected_material_config",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            [
                {"element": "Si", "coordinate_2d": [0.5, 0.5], "distance_z": 2.0},
                {"element": "C", "coordinate_2d": [0.25, 0.25], "distance_z": 1.5},
                {"element": "N", "coordinate_2d": [0.75, 0.75], "distance_z": 1.0},
            ],
            AdatomPlacementMethodEnum.EXACT_COORDINATE.value,
            SLAB_Si_3_ADATOMS,
        ),
    ],
)
def test_create_multiple_adatom_defects(crystal_config, adatom_dicts, placement_method, expected_material_config):
    slab = MaterialWithBuildMetadata.create(crystal_config)

    defect_models = []
    for adatom_data in adatom_dicts:
        defect_dict = AdatomDefectDict(**adatom_data)
        defect_models.append(defect_dict)

    defects = create_multiple_adatom_defects(slab, defect_models, placement_method)
    defects.metadata.build = []
    assert_two_entities_deep_almost_equal(defects, expected_material_config)
