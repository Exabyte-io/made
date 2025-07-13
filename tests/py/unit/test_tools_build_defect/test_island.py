from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.slab.helpers import create_island_defect
from mat3ra.made.tools.build.slab.helpers import create_slab
from mat3ra.made.tools.utils import coordinate as CoordinateCondition
from unit.fixtures.bulk import BULK_Si_CONVENTIONAL


# @pytest.mark.parametrize(
#     "crystal_config, condition_class, condition_params, num_added_layers, num_atoms_in_island, expected_last_element",
#     [
#         (
#             SI_CONVENTIONAL_SLAB_001,
#             CoordinateCondition.CylinderCoordinateCondition,
#             {"center_position": [0.5, 0.5], "radius": 0.55, "min_z": 0, "max_z": 1},
#             1,
#             1,
#             "Si",
#         )
#     ],
# )
def test_create_island_defect_new_pattern():
    slab = create_slab(
        Material.create(BULK_Si_CONVENTIONAL),
        xy_supercell_matrix=[[3, 0], [0, 3]],
    )
    condition = CoordinateCondition.CylinderCoordinateCondition(
        center_position=[0.5, 0.5], radius=0.25, min_z=0, max_z=1
    )

    defect = create_island_defect(
        slab=slab,
        condition=condition,
        use_cartesian_coordinates=False,
        number_of_added_layers=1,
    )

    assert defect is not None
