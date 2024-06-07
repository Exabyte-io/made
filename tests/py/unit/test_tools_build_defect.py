from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import PointDefectConfiguration, create_defect


def test_create_defect():
    # vacancy in place of 0 element
    configuration = PointDefectConfiguration()
    material = Material.create(Material.default_config)
    defect = create_defect(material, configuration)

    assert len(defect.basis["elements"]) == 1
