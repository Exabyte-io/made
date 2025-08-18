import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defective_structures.two_dimensional.island.helpers import create_defect_island
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.helpers import create_slab
from mat3ra.made.tools.entities import coordinate as CoordinateCondition
from unit.fixtures.bulk import BULK_Si_CONVENTIONAL
from unit.fixtures.island import ISLAND_SLAB_Si_001_CYLINDER_CONDITION_1_LAYER
from unit.utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "slab_parameters, condition_class, condition_params, num_added_layers, expected_config",
    [
        (
            {"crystal": BULK_Si_CONVENTIONAL, "number_of_layers": 2, "xy_supercell_matrix": [[3, 0], [0, 3]]},
            CoordinateCondition.CylinderCoordinateCondition,
            {"center_position": [0.5, 0.5], "radius": 0.25, "min_z": 0, "max_z": 1},
            1,
            ISLAND_SLAB_Si_001_CYLINDER_CONDITION_1_LAYER,
        )
    ],
)
def test_create_island_defect(slab_parameters, condition_class, condition_params, num_added_layers, expected_config):
    slab = create_slab(
        Material.create(slab_parameters["crystal"]),
        number_of_layers=slab_parameters["number_of_layers"],
        xy_supercell_matrix=slab_parameters["xy_supercell_matrix"],
        termination_top_formula=None,
        termination_bottom_formula=None,
    )
    condition = condition_class(
        center_position=condition_params["center_position"],
        radius=condition_params["radius"],
        min_z=condition_params["min_z"],
        max_z=condition_params["max_z"],
    )

    defect = create_defect_island(
        slab=slab,
        condition=condition,
        use_cartesian_coordinates=False,
        number_of_added_layers=num_added_layers,
    )
    defect.metadata.build = []
    assert_two_entities_deep_almost_equal(defect, expected_config)
