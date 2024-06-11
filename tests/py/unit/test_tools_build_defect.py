from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import (
    InterstitialConfiguration,
    SubstitutionConfiguration,
    VacancyConfiguration,
    create_defect,
)


def test_create_vacancy():
    # vacancy in place of 0 element
    configuration = VacancyConfiguration()
    material = Material.create(Material.default_config)
    defect = create_defect(material, configuration)

    assert len(defect.basis["elements"]) == 1


def test_create_substitution():
    # Substitution of Ge in place of Si at 0 element
    configuration = SubstitutionConfiguration(element="Ge")
    material = Material.create(Material.default_config)
    defect = create_defect(material, configuration)

    assert defect.basis["elements"] == [{"id": 0, "value": "Ge"}, {"id": 1, "value": "Si0+"}]
    assert defect.basis["coordinates"][0] == {"id": 0, "value": [0.0, 0.0, 0.0]}


def test_create_interstitial():
    # Interstitial Ge at 0.5, 0.5, 0.5
    configuration = InterstitialConfiguration(element="Ge", position=[0.5, 0.5, 0.5])
    material = Material.create(Material.default_config)
    defect = create_defect(material, configuration)

    assert defect.basis["elements"] == [
        {"id": 0, "value": "Ge3+"},
        {"id": 1, "value": "Si0+"},
        {"id": 2, "value": "Si0+"},
    ]
