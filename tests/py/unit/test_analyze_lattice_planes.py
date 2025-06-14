import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.build.slab.entities import Termination
from mat3ra.made.tools.convert import from_pymatgen
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.generated.fixtures import HfO2_BULK_MATERIAL, SrTiO3_BULK_MATERIAL

MILLER_INDICES = [(0, 0, 1), (1, 1, 0), (0, 1, 1)]
SrTiO3_EXPECTED_TERMINATION_FORMULAS = {
    (0, 0, 1): ["SrO", "TiO2"],
    (1, 1, 0): ["O2", "SrTiO"],
    (0, 1, 1): ["O2", "SrTiO"],
}
SrTiO3_EXPECTED_TERMINATIONS_WITH_VACUUM = {
    (0, 0, 1): [
        Termination.from_string("SrO_Pmmm_2"),
        Termination.from_string("TiO2_Pmmm_2"),
    ],
    (1, 1, 0): [
        Termination.from_string("O2_Pmmm_2"),
        Termination.from_string("SrTiO_Pmmm_3"),
    ],
    (0, 1, 1): [
        Termination.from_string("O2_Pmmm_2"),
        Termination.from_string("SrTiO_Pmmm_3"),
    ],
}

SrTiO3_EXPECTED_TERMINATIONS_WITHOUT_VACUUM = {
    (0, 0, 1): [
        None,
        Termination.from_string("SrO_P4/mmm_2"),
    ],
    (1, 1, 0): [
        None,
        Termination.from_string("SrTiO_P4/mmm_3"),
    ],
    (0, 1, 1): [
        None,
        Termination.from_string("SrTiO_P4/mmm_3"),
    ],
}

# 3 slabs
SrTiO3_SLAB_WITH_VACUUM_NAMES = ["Sr3 Ti3 O7", "Sr2 Ti2 O8"]
# Only one possibility
SrTiO3_SLAB_WITHOUT_VACUUM_NAMES = [
    "Sr1 Ti1 O3",
]
# shifts for 3 slabs
SrTiO3_SHIFTS_WITH_VACUUM = [0.25, 0.25]
# Only one possibility
SrTiO3_SHIFTS_WITHOUT_VACUUM = [0.0, 0.25]

# shifts for 3 slabs
HfO2_SHIFTS_WITH_VACUUM = [0.0, 0.25, 0.25]
HfO2_SHIFTS_WITHOUT_VACUUM = [0.0, 0.0]
# shifts for 3 slabs
HfO2_SLAB_WITH_VACUUM_NAMES = [
    "Hf8 O22",
    "Hf8 O20",
    "Hf8 O18",
    "Hf12 O24",
    "Hf12 O22",
    "Hf12 O20",
    "Hf12 O18",
    "Hf16 O24",
]
HfO2_SLAB_WITHOUT_VACUUM_NAMES = ["Hf4 O8"]


@pytest.mark.parametrize(
    "material, miller_indices, expected_slabs_with_vacuum_names",
    [
        (SrTiO3_BULK_MATERIAL, (1, 1, 0), SrTiO3_SLAB_WITH_VACUUM_NAMES),
        (HfO2_BULK_MATERIAL, (1, 1, 1), HfO2_SLAB_WITH_VACUUM_NAMES),
    ],
)
def test_all_planes_as_pymatgen_slabs_with_vacuum(material, miller_indices, expected_slabs_with_vacuum_names):
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)

    slabs_with_vacuum = [
        Material.create(from_pymatgen(slab)) for slab in analyzer.all_planes_as_pymatgen_slabs_with_vacuum
    ]

    assert len(slabs_with_vacuum) == len(expected_slabs_with_vacuum_names)
    assert all(slab.name in expected_slabs_with_vacuum_names for slab in slabs_with_vacuum)


@pytest.mark.parametrize(
    "material, miller_indices, expected_slabs_without_vacuum_names",
    [
        (SrTiO3_BULK_MATERIAL, (1, 1, 0), SrTiO3_SLAB_WITHOUT_VACUUM_NAMES),
        (HfO2_BULK_MATERIAL, (1, 1, 1), HfO2_SLAB_WITHOUT_VACUUM_NAMES),
    ],
)
def test_all_planes_as_pymatgen_slabs_without_vacuum(material, miller_indices, expected_slabs_without_vacuum_names):
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)

    slabs_without_vacuum = [
        Material.create(from_pymatgen(slab)) for slab in analyzer.all_planes_as_pymatgen_slabs_without_vacuum
    ]

    assert len(slabs_without_vacuum) == len(expected_slabs_without_vacuum_names)
    assert all(slab.name in expected_slabs_without_vacuum_names for slab in slabs_without_vacuum)


@pytest.mark.parametrize(
    "material, miller_indices, expected_shifts_with_vacuum",
    [
        (SrTiO3_BULK_MATERIAL, (1, 1, 0), SrTiO3_SHIFTS_WITH_VACUUM),
        (HfO2_BULK_MATERIAL, (1, 1, 1), HfO2_SHIFTS_WITH_VACUUM),
    ],
)
def test_termination_holders_shifts_with_vacuum(material, miller_indices, expected_shifts_with_vacuum):
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    shifts_with_vacuum = [h.shift_with_vacuum for h in analyzer.termination_holders]
    assertion_utils.assert_deep_almost_equal(shifts_with_vacuum, expected_shifts_with_vacuum)


@pytest.mark.parametrize(
    "material, miller_indices, expected_shifts_without_vacuum",
    [
        (SrTiO3_BULK_MATERIAL, (1, 1, 0), SrTiO3_SHIFTS_WITHOUT_VACUUM),
        (HfO2_BULK_MATERIAL, (1, 1, 1), HfO2_SHIFTS_WITHOUT_VACUUM),
    ],
)
def test_termination_holders_shifts_without_vacuum(material, miller_indices, expected_shifts_without_vacuum):
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    shifts_without_vacuum = [
        h.shift_without_vacuum for h in analyzer.termination_holders if h.shift_without_vacuum is not None
    ]
    assertion_utils.assert_deep_almost_equal(shifts_without_vacuum, expected_shifts_without_vacuum)


@pytest.mark.parametrize(
    "miller_indices",
    [(0, 0, 1), (1, 1, 0), (0, 1, 1)],
)
def test_terminations_with_vacuum(miller_indices):
    material = SrTiO3_BULK_MATERIAL
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    terminations_with_vacuum = analyzer.terminations_with_vacuum
    terminations_with_vacuum_str_list = [str(t) for t in terminations_with_vacuum]
    expected_with_vacuum = [str(t) for t in SrTiO3_EXPECTED_TERMINATIONS_WITH_VACUUM[miller_indices]]
    assert terminations_with_vacuum_str_list == expected_with_vacuum


@pytest.mark.parametrize(
    "miller_indices",
    [(0, 0, 1), (1, 1, 0), (0, 1, 1)],
)
def test_terminations_without_vacuum(miller_indices):
    material = SrTiO3_BULK_MATERIAL
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    terminations_without_vacuum = analyzer.terminations_without_vacuum
    terminations_without_vacuum_str_list = [str(t) for t in terminations_without_vacuum]
    expected_without_vacuum = [str(t) for t in SrTiO3_EXPECTED_TERMINATIONS_WITHOUT_VACUUM[miller_indices]]
    assert terminations_without_vacuum_str_list == expected_without_vacuum


@pytest.mark.parametrize("miller_indices", MILLER_INDICES)
def test_terminations(miller_indices):
    material = SrTiO3_BULK_MATERIAL
    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    terminations = analyzer.terminations
    actual_formulas = sorted([t.formula for t in terminations])
    expected_formulas = sorted(SrTiO3_EXPECTED_TERMINATION_FORMULAS[miller_indices])
    assert actual_formulas == expected_formulas


def test_110_and_011_terminations_for_srtio3_are_the_same():
    analyzer_110 = CrystalLatticePlanesMaterialAnalyzer(material=SrTiO3_BULK_MATERIAL, miller_indices=(1, 1, 0))
    analyzer_011 = CrystalLatticePlanesMaterialAnalyzer(material=SrTiO3_BULK_MATERIAL, miller_indices=(0, 1, 1))
    formulas_110 = sorted([t.formula for t in analyzer_110.terminations])
    formulas_011 = sorted([t.formula for t in analyzer_011.terminations])
    assert formulas_110 == formulas_011
