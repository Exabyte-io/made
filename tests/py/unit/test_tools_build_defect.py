from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import (
    AdatomSlabPointDefectConfiguration,
    CrystalSiteAdatomSlabDefectBuilder,
    EquidistantAdatomSlabDefectBuilder,
    PointDefectBuilderParameters,
    PointDefectConfiguration,
    create_defect,
    create_slab_defect,
)
from mat3ra.made.tools.build.defect.builders import IslandSlabDefectBuilder
from mat3ra.made.tools.build.defect.configuration import IslandSlabDefectConfiguration
from mat3ra.made.tools.build.slab import SlabConfiguration, create_slab, get_terminations
from mat3ra.made.tools.utils import CoordinateCondition
from mat3ra.utils import assertion as assertion_utils

clean_material = Material.create(Material.default_config)

slab_config = SlabConfiguration(
    bulk=clean_material,
    miller_indices=(1, 1, 1),
    thickness=4,
    vacuum=6,
    xy_supercell_matrix=[[1, 0], [0, 1]],
    use_orthogonal_z=True,
)
t = get_terminations(slab_config)[0]
slab = create_slab(slab_config, t)


def test_create_vacancy():
    # vacancy in place of 0 element
    configuration = PointDefectConfiguration(crystal=clean_material, defect_type="vacancy", site_id=0)
    defect = create_defect(configuration)

    assert len(defect.basis.elements.values) == 1


def test_create_substitution():
    # Substitution of Ge in place of Si at default site_id=0
    configuration = PointDefectConfiguration(crystal=clean_material, defect_type="substitution", chemical_element="Ge")
    defect = create_defect(configuration)

    assert defect.basis.elements.to_dict() == [{"id": 0, "value": "Ge"}, {"id": 1, "value": "Si"}]
    assert defect.basis.coordinates.to_dict()[0] == {"id": 0, "value": [0.0, 0.0, 0.0]}


def test_create_interstitial():
    # Interstitial Ge at 0.5, 0.5, 0.5 position
    configuration = PointDefectConfiguration(
        crystal=clean_material, defect_type="interstitial", chemical_element="Ge", position=[0.5, 0.5, 0.5]
    )
    defect = create_defect(configuration)

    assert defect.basis.elements.to_dict() == [
        {"id": 0, "value": "Ge"},
        {"id": 1, "value": "Si"},
        {"id": 2, "value": "Si"},
    ]


def test_create_defect_from_site_id():
    # Substitution of Ge in place of Si at site_id=1
    defect_configuration = PointDefectConfiguration.from_site_id(
        crystal=clean_material, defect_type="substitution", chemical_element="Ge", site_id=1
    )
    defect_builder_parameters = PointDefectBuilderParameters(center_defect=False)
    material_with_defect = create_defect(
        builder_parameters=defect_builder_parameters, configuration=defect_configuration
    )

    assert material_with_defect.basis.elements.to_dict() == [
        {"id": 0, "value": "Si"},
        {"id": 1, "value": "Ge"},
    ]


def test_create_adatom():
    # Adatom of Si at 0.5, 0.5 position
    configuration = AdatomSlabPointDefectConfiguration(
        crystal=slab, position_on_surface=[0.5, 0.5], distance_z=2, chemical_element="Si"
    )
    defect = create_slab_defect(configuration=configuration, builder=None)

    assert defect.basis.elements.values[-1] == "Si"
    assertion_utils.assert_deep_almost_equal([0.5, 0.5, 0.450843412], defect.basis.coordinates.values[-1])


def test_create_adatom_equidistant():
    # Adatom of Si at approximate 0.5, 0.5 position
    configuration = AdatomSlabPointDefectConfiguration(
        crystal=slab, position_on_surface=[0.5, 0.5], distance_z=2, chemical_element="Si"
    )
    defect = create_slab_defect(configuration=configuration, builder=EquidistantAdatomSlabDefectBuilder())

    assert defect.basis.elements.values[-1] == "Si"
    # We expect adatom to shift from provided position
    assertion_utils.assert_deep_almost_equal(
        [0.583333334, 0.458333333, 0.450843412], defect.basis.coordinates.values[-1]
    )


def test_create_crystal_site_adatom():
    # Adatom of Si (autodetect) at approximate 0.5, 0.5 position
    configuration = AdatomSlabPointDefectConfiguration(
        crystal=slab, position_on_surface=[0.5, 0.5], distance_z=2, chemical_element=None
    )
    builder = CrystalSiteAdatomSlabDefectBuilder()
    defect = create_slab_defect(configuration=configuration, builder=builder)

    assert defect.basis.elements.values[-1] == "Si"
    assertion_utils.assert_deep_almost_equal(
        [0.083333333, 0.458333333, 0.352272727], defect.basis.coordinates.values[-1]
    )


def test_create_island():
    ni_config = {
        "lattice": {
            "a": 3.475146,
            "b": 3.475146,
            "c": 3.475146,
            "alpha": 90,
            "beta": 90,
            "gamma": 120,
            "units": {"length": "angstrom", "angle": "degree"},
            "type": "HEX",
            "vectors": {
                "a": [3.475146, 0, 0],
                "b": [-1.737573, 3.009565, 0],
                "c": [0, 0, 3.475146],
                "alat": 1,
                "units": "angstrom",
            },
        },
        "basis": {
            "elements": [
                {"id": 0, "value": "Ni"},
                {"id": 1, "value": "Ni"},
                {"id": 2, "value": "Ni"},
                {"id": 3, "value": "Ni"},
            ],
            "coordinates": [
                {"id": 0, "value": [0, 0, 0]},
                {"id": 1, "value": [0, 0.5, 0.5]},
                {"id": 2, "value": [0.5, 0, 0.5]},
                {"id": 3, "value": [0.5, 0.5, 0]},
            ],
            "units": "crystal",
            "cell": [
                [3.475146, 0, 2.127913212734864e-16],
                [-1.7375729999999994, 3.0095647178598774, 2.127913212734864e-16],
                [0, 0, 3.475146],
            ],
            "constraints": [],
        },
        "name": "Ni HEX",
        "isNonPeriodic": False,
    }

    clean_material = Material.create(ni_config)

    print("clean_material", clean_material.basis.coordinates.values)

    slab_config = SlabConfiguration(
        bulk=clean_material,
        miller_indices=(1, 1, 1),
        thickness=4,
        vacuum=6,
        xy_supercell_matrix=[[1, 0], [0, 1]],
        use_orthogonal_z=True,
        make_primitive=False,
    )
    print(get_terminations(slab_config))
    t = get_terminations(slab_config)[0]
    slab = create_slab(slab_config, t)
    print(slab.to_json())
    print("slab", slab.basis.coordinates.values)

    condition = CoordinateCondition.cylinder(center_position=[0.5, 0.5], radius=0.25, min_z=0, max_z=1)
    island_config = IslandSlabDefectConfiguration(
        crystal=slab,
        defect_type="island",
        condition=condition,
        thickness=1,
    )

    defect = create_slab_defect(configuration=island_config, builder=IslandSlabDefectBuilder())
    print("slab elements", len(slab.basis.elements.values), slab.basis.elements.values)
    print("slab coords", len(slab.basis.coordinates.values), slab.basis.coordinates.values)
    print("defect elements", len(defect.basis.elements.values), defect.basis.elements.values)
    print("defect coords", len(defect.basis.coordinates.values), defect.basis.coordinates.values)
    # Only one atom is in the island for this configuration
    assert len(defect.basis.elements.values) == len(slab.basis.elements.values) + 10
    assert defect.basis.elements.values[-1] == "Si"
