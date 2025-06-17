import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import create_slab_defect
from mat3ra.made.tools.build.defect.builders import IslandSlabDefectBuilder
from mat3ra.made.tools.build.defect.configuration import IslandSlabDefectConfiguration
from mat3ra.made.tools.utils import coordinate as CoordinateCondition
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001


@pytest.mark.skip(reason="we'll fix before epic-7623 is merged")
@pytest.mark.parametrize(
    "crystal_config, condition_params, defect_type, num_added_layers, num_atoms_in_island, expected_last_element",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
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
