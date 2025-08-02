import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.bulk import BULK_SrTiO3

SLAB_SrTiO3_EXPECTED_TERMINATIONS_WITH_VACUUM = [
    [
        (1, 1, 0),
        ["SrTiO_Pmmm_3", "O2_Pmmm_2"],
    ]
]

SLAB_SrTiO3_EXPECTED_TERMINATIONS_WITHOUT_VACUUM = [
    [
        (1, 1, 0),
        [
            "None",
            "SrTiO_P4/mmm_3",
        ],
    ]
]


@pytest.mark.parametrize(
    "material_config, miller_indices, number_of_slabs",
    [
        (
            BULK_SrTiO3,
            (1, 1, 0),
            2,
        ),
    ],
)
def test_all_planes_as_pymatgen_slabs_with_vacuum(material_config, miller_indices, number_of_slabs):
    material = Material.create(material_config)
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    slabs_with_vacuum = analyzer.all_planes_as_pymatgen_slabs_with_vacuum
    assert len(slabs_with_vacuum) == number_of_slabs


@pytest.mark.parametrize(
    "material_config, miller_indices, number_of_slabs",
    [
        (BULK_SrTiO3, (1, 1, 0), 1),
    ],
)
def test_all_planes_as_pymatgen_slabs_without_vacuum(material_config, miller_indices, number_of_slabs):
    material = Material.create(material_config)
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    slabs_without_vacuum = analyzer.all_planes_as_pymatgen_slabs_without_vacuum
    assert len(slabs_without_vacuum) == number_of_slabs


@pytest.mark.parametrize(
    "material_config, miller_indices, expected_shifts_with_vacuum",
    [
        (BULK_SrTiO3, (1, 1, 0), [0.25, 0.25]),
    ],
)
def test_termination_holders_shifts_with_vacuum(material_config, miller_indices, expected_shifts_with_vacuum):
    material = Material.create(material_config)
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    shifts_with_vacuum = [h.shift_with_vacuum for h in analyzer.termination_holders]
    assertion_utils.assert_deep_almost_equal(shifts_with_vacuum, expected_shifts_with_vacuum)


@pytest.mark.parametrize(
    "material_config, miller_indices, expected_shifts_without_vacuum",
    [
        (BULK_SrTiO3, (1, 1, 0), [0.0, 0.25]),
    ],
)
def test_termination_holders_shifts_without_vacuum(material_config, miller_indices, expected_shifts_without_vacuum):
    material = Material.create(material_config)
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    shifts_without_vacuum = [
        h.shift_without_vacuum for h in analyzer.termination_holders if h.shift_without_vacuum is not None
    ]
    assertion_utils.assert_deep_almost_equal(shifts_without_vacuum, expected_shifts_without_vacuum)


@pytest.mark.parametrize(
    "material_config, miller_indices, expected_terminations_with_vacuum_str",
    [
        (
            BULK_SrTiO3,
            (1, 1, 0),
            SLAB_SrTiO3_EXPECTED_TERMINATIONS_WITH_VACUUM[0][1],
        )
    ],
)
def test_terminations_with_vacuum(material_config, miller_indices, expected_terminations_with_vacuum_str):
    material = Material.create(material_config)
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    terminations_with_vacuum = analyzer.terminations_with_vacuum
    terminations_with_vacuum_str_list = [str(t) for t in terminations_with_vacuum]
    assert sorted(terminations_with_vacuum_str_list) == sorted(expected_terminations_with_vacuum_str)


@pytest.mark.parametrize(
    "material_config, miller_indices, expected_terminations_without_vacuum_str",
    [
        (
            BULK_SrTiO3,
            (1, 1, 0),
            SLAB_SrTiO3_EXPECTED_TERMINATIONS_WITHOUT_VACUUM[0][1],
        ),
    ],
)
def test_terminations_without_vacuum(material_config, miller_indices, expected_terminations_without_vacuum_str):
    material = Material.create(material_config)
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    terminations_without_vacuum = analyzer.terminations_without_vacuum
    terminations_without_vacuum_str_list = [str(t) for t in terminations_without_vacuum]
    all_materials = analyzer.get_materials_for_all_terminations_without_vacuum()
    assert sorted(terminations_without_vacuum_str_list) == sorted(expected_terminations_without_vacuum_str)
    assert all_materials
    assert len(all_materials) == len([t for t in analyzer.terminations_without_vacuum if t is not None])


@pytest.mark.parametrize(
    "material_config, miller_indices, expected_formulas",
    [(BULK_SrTiO3, (1, 1, 0), ["SrTiO", "O2"])],
)
def test_terminations(material_config, miller_indices, expected_formulas):
    material = Material.create(material_config)
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    terminations = analyzer.terminations
    actual_formulas = sorted([t.formula for t in terminations])
    assert sorted(actual_formulas) == sorted(expected_formulas)
