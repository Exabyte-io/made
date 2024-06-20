from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import PointDefectBuilderParameters, PointDefectConfiguration, create_defect

clean_material = Material.create(Material.default_config)


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
