import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.slab.configuration import VoidSite
from mat3ra.made.tools.build.defect.slab.helpers import create_island_defect
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.utils import coordinate as CoordinateCondition

from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001


def test_create_island_defect_new_pattern():
    slab = Material.create(SI_CONVENTIONAL_SLAB_001)
    crystal = create_supercell(slab, scaling_factor=[3, 3, 1])
    condition = CoordinateCondition.CylinderCoordinateCondition(
        center_position=[0.5, 0.5], radius=0.15, min_z=0, max_z=1
    )

    defect = create_island_defect(
        slab=crystal,
        condition=condition,
        number_of_added_layers=1,
    )

    assert defect is not None


@pytest.mark.parametrize(
    "condition_class, condition_params",
    [
        (
            CoordinateCondition.CylinderCoordinateCondition,
            {"center_position": [0.5, 0.5], "radius": 0.15, "min_z": 0, "max_z": 1},
            "island",
            1,
            1,
            "Si",
        )
    ],
)
def test_create_island(
    crystal_config, condition_params, defect_type, num_added_layers, num_atoms_in_island, expected_last_element
):
    # TODO: use TiN
    crystal = Material.create(crystal_config)
    condition = CoordinateCondition.CylinderCoordinateCondition(**condition_params)
    island_config = IslandSlabDefectConfiguration(
        crystal=crystal,
        defect_type=defect_type,
        condition=condition,
        number_of_added_layers=num_added_layers,
    )

    defect = create_slab_defect(configuration=island_config, builder=IslandSlabDefectBuilder())

    # 1 atom in the island were added for this configuration with 001 slab orientation
    assert len(defect.basis.elements.values) == len(crystal.basis.elements.values) + num_atoms_in_island
    assert defect.basis.elements.values[-1] == expected_last_element
