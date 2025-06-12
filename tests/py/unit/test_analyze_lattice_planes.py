from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.build.slab.entities import Termination, TerminationHolder
from mat3ra.made.tools.convert import from_pymatgen
from unit.fixtures.generated.fixtures import SrTiO3_BULK_MATERIAL, HfO2_BULK_MATERIAL

from .utils import assert_two_entities_deep_almost_equal

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
        Termination.from_string("SrO_Pmmm_2"),
        Termination.from_string("TiO2_Pmmm_2"),
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

# TODO: adjust to match real results (13 terminations)
# shifts for 3 slabs
HfO2_SHIFTS_WITH_VACUUM = [0.0, 0.25, 0.25]
HfO2_SHIFTS_WITHOUT_VACUUM = [0.0, 0.0]
# shifts for 3 slabs
HfO2_SLAB_WITH_VACUUM_NAMES = ["Hf1 O2", "Hf2 O4", "Hf3 O6"]
HfO2_SLAB_WITHOUT_VACUUM_NAMES = ["Hf1 O2", "Hf1 O2"]


def create_analyzer(material, miller_indices):
    return CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)


def test_crystal_lattice_planes_analyzer_slabs_110():
    analyzer = create_analyzer(SrTiO3_BULK_MATERIAL, (1, 1, 0))

    slabs_without_vacuum = [
        Material.create(from_pymatgen(slab)) for slab in analyzer.all_planes_as_pymatgen_slabs_without_vacuum
    ]
    slabs_with_vacuum = [
        Material.create(from_pymatgen(slab)) for slab in analyzer.all_planes_as_pymatgen_slabs_with_vacuum
    ]

    assert len(slabs_with_vacuum) == len(SrTiO3_SLAB_WITH_VACUUM_NAMES)
    assert len(slabs_without_vacuum) == len(SrTiO3_SLAB_WITHOUT_VACUUM_NAMES)
    assert all(slab.name in SrTiO3_SLAB_WITH_VACUUM_NAMES for slab in slabs_with_vacuum)
    assert all(slab.name in SrTiO3_SLAB_WITHOUT_VACUUM_NAMES for slab in slabs_without_vacuum)


def test_crystal_lattice_planes_analyzer_shifts_011():
    analyzer = create_analyzer(SrTiO3_BULK_MATERIAL, (1, 1, 0))
    shifts_with_vacuum = [h.shift_with_vacuum for h in analyzer.termination_holders]
    shifts_without_vacuum = [
        h.shift_without_vacuum for h in analyzer.termination_holders if h.shift_without_vacuum is not None
    ]

    assert shifts_with_vacuum == SrTiO3_SHIFTS_WITH_VACUUM
    assert shifts_without_vacuum == SrTiO3_SHIFTS_WITHOUT_VACUUM


def test_crystal_lattice_planes_analyzer_terminations_011():
    analyzer = create_analyzer(SrTiO3_BULK_MATERIAL, (1, 1, 0))
    terminations_with_vacuum = [h.termination_with_vacuum for h in analyzer.termination_holders]
    terminations_without_vacuum = [h.termination_without_vacuum for h in analyzer.termination_holders]

    terminations_with_vacuum_str_list = [str(t) for t in terminations_with_vacuum]
    terminations_without_vacuum_str_list = [str(t) for t in terminations_without_vacuum]
    assert terminations_with_vacuum_str_list == [str(t) for t in SrTiO3_EXPECTED_TERMINATIONS_WITH_VACUUM[(1, 1, 0)]]
    assert terminations_without_vacuum_str_list == [
        str(t) for t in SrTiO3_EXPECTED_TERMINATIONS_WITHOUT_VACUUM[(1, 1, 0)]
    ]


def test_crystal_lattice_planes_analyzer_slabs_111():
    analyzer = create_analyzer(HfO2_BULK_MATERIAL, (1, 1, 1))

    slabs_without_vacuum = [
        Material.create(from_pymatgen(slab)) for slab in analyzer.all_planes_as_pymatgen_slabs_without_vacuum
    ]
    slabs_with_vacuum = [
        Material.create(from_pymatgen(slab)) for slab in analyzer.all_planes_as_pymatgen_slabs_with_vacuum
    ]

    assert len(slabs_with_vacuum) == len(HfO2_SLAB_WITH_VACUUM_NAMES)
    assert len(slabs_without_vacuum) == len(HfO2_SLAB_WITHOUT_VACUUM_NAMES)
    assert all(slab.name in HfO2_SLAB_WITH_VACUUM_NAMES for slab in slabs_with_vacuum)
    assert all(slab.name in HfO2_SLAB_WITHOUT_VACUUM_NAMES for slab in slabs_without_vacuum)


def test_terminations():
    for miller_indices in MILLER_INDICES:
        analyzer = create_analyzer(SrTiO3_BULK_MATERIAL, miller_indices)
        terminations = analyzer.terminations
        actual_formulas = sorted([t.formula for t in terminations])
        expected_formulas = sorted(SrTiO3_EXPECTED_TERMINATION_FORMULAS[miller_indices])
        assert actual_formulas == expected_formulas


def test_110_and_011_terminations_for_srtio3_are_the_same():
    analyzer_110 = create_analyzer(SrTiO3_BULK_MATERIAL, (1, 1, 0))
    analyzer_011 = create_analyzer(SrTiO3_BULK_MATERIAL, (0, 1, 1))
    formulas_110 = sorted([t.formula for t in analyzer_110.terminations])
    formulas_011 = sorted([t.formula for t in analyzer_011.terminations])
    assert formulas_110 == formulas_011
