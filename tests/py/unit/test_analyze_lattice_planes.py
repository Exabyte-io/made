from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice import LatticeMaterialAnalyzer
from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.build.slab.entities import TerminationHolder, Termination
from mat3ra.made.tools.convert import from_pymatgen
from unit.fixtures.generated.fixtures import SrTiO3_BULK_MATERIAL
from .utils import assert_two_entities_deep_almost_equal

MATERIAL = SrTiO3_BULK_MATERIAL
MILLER_INDICES = [(0, 0, 1), (1, 1, 0), (0, 1, 1)]
EXPECTED_TERMINATIONS = {
    (0, 0, 1): ["SrO", "TiO2"],
    (1, 1, 0): ["O2", "SrTiO"],
    (0, 1, 1): ["O2", "SrTiO"],
}


def create_analyzer(material, miller_indices):
    return CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)


def test_crystal_lattice_planes_analyzer_011():
    conventional_material = LatticeMaterialAnalyzer(material=MATERIAL).material_with_conventional_lattice
    analyzer = create_analyzer(conventional_material, (1, 1, 0))

    slabs_without_vacuum = [
        Material.create(from_pymatgen(slab)) for slab in analyzer.all_planes_as_pymatgen_slabs_without_vacuum
    ]
    slabs_with_vacuum = [
        Material.create(from_pymatgen(slab)) for slab in analyzer.all_planes_as_pymatgen_slabs_with_vacuum
    ]
    shifts_without_vacuum = [slab.shift for slab in analyzer.all_planes_as_pymatgen_slabs_with_vacuum]
    shifts_with_vacuum = [slab.shift for slab in analyzer.all_planes_as_pymatgen_slabs_without_vacuum]

    # 3 slabs
    SLAB_WITH_VACUUM_NAMES = ["Sr3 Ti3 O7", "Sr2 Ti2 O8", "Sr3 Ti3 O9"]

    # Only one possibility
    SLAB_WITHOUT_VACUUM_NAMES = [
        "Sr1 Ti1 O3",
    ]

    # shifts for 3 slabs
    SHIFTS_WITH_VACUUM = [0.0, 0.25, 0.25]

    # Only one possibility
    SHIFTS_WITHOUT_VACUUM = [0.0]

    assert all(slab.name in SLAB_WITH_VACUUM_NAMES for slab in slabs_with_vacuum)
    assert all(slab.name in SLAB_WITHOUT_VACUUM_NAMES for slab in slabs_without_vacuum)
    assert shifts_with_vacuum == SHIFTS_WITH_VACUUM
    assert shifts_without_vacuum == SHIFTS_WITHOUT_VACUUM


def test_termination_holders_for_011():
    MILLER_INDICES = (0, 1, 1)
    analyzer = create_analyzer(MATERIAL, (0, 1, 1))
    termination_holders = analyzer.termination_holders

    expected_data = TerminationHolder(
        termination_with_vacuum=Termination.from_string("SrTiO"),
        termination_without_vacuum=Termination.from_string(EXPECTED_TERMINATIONS["TiO2"]),
        shift_with_vacuum=Termination.from_string(EXPECTED_TERMINATIONS["SrTiO"]),
        shift_without_vacuum=Termination.from_string(EXPECTED_TERMINATIONS["O2"]),
    )

    assert_two_entities_deep_almost_equal(termination_holders, expected_data)


def test_terminations():
    for miller_indices in MILLER_INDICES:
        analyzer = create_analyzer(MATERIAL, miller_indices)
        terminations = analyzer.terminations
        actual_formulas = sorted([t.formula for t in terminations])
        expected_formulas = sorted(EXPECTED_TERMINATIONS[miller_indices])
        assert actual_formulas == expected_formulas, f"Mismatch for MI {miller_indices}"


def test_110_and_011_terminations_are_the_same():
    analyzer_110 = create_analyzer(MATERIAL, (1, 1, 0))
    analyzer_011 = create_analyzer(MATERIAL, (0, 1, 1))
    formulas_110 = sorted([t.formula for t in analyzer_110.terminations])
    formulas_011 = sorted([t.formula for t in analyzer_011.terminations])
    assert formulas_110 == formulas_011


def test_termination_holders():
    EXPECTED_HOLDERS_DATA = {
        (0, 0, 1): [
            {"termination_with_vacuum": "SrO", "termination_without_vacuum": "SrO", "shift": 0.5},
            {"termination_with_vacuum": "TiO2", "termination_without_vacuum": "TiO2", "shift": 0.0},
        ],
        (1, 1, 0): [
            {"termination_with_vacuum": "SrTiO", "termination_without_vacuum": "SrTiO", "shift": 0.0},
            {"termination_with_vacuum": "O2", "termination_without_vacuum": "O2", "shift": 0.25},
        ],
        (0, 1, 1): [
            {"termination_with_vacuum": "SrTiO", "termination_without_vacuum": "SrTiO", "shift": 0.0},
            {"termination_with_vacuum": "O2", "termination_without_vacuum": "O2", "shift": 0.25},
        ],
    }

    for miller_indices, expected_data in EXPECTED_HOLDERS_DATA.items():
        analyzer = create_analyzer(MATERIAL, miller_indices)
        termination_holders = analyzer.termination_holders

        actual_holders_as_dict = [
            {
                "termination_with_vacuum": h.termination_with_vacuum.formula,
                "termination_without_vacuum": h.termination_without_vacuum.formula,
                "shift": h.shift,
            }
            for h in termination_holders
        ]

        # Sorting by a consistent key to ensure order doesn't matter
        sorted_actual = sorted(actual_holders_as_dict, key=lambda x: x["termination_with_vacuum"])
        sorted_expected = sorted(expected_data, key=lambda x: x["termination_with_vacuum"])

        assert len(sorted_actual) == len(sorted_expected)
        for actual, expected_holder in zip(sorted_actual, sorted_expected):
            assert actual["termination_with_vacuum"] == expected_holder["termination_with_vacuum"]
            assert actual["termination_without_vacuum"] == expected_holder["termination_without_vacuum"]
            assert abs(actual["shift"] - expected_holder["shift"]) < 1e-4, f"Shift mismatch for {miller_indices}"
