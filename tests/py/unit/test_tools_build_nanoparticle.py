import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.pristine_structures.zero_dimensional.nanoparticle import (
    NanoparticleShapesEnum,
    create_nanoparticle_by_shape,
    create_nanoparticle_from_material,
)
from mat3ra.made.tools.entities.coordinate import SphereCoordinateCondition
from unit.fixtures.bulk import BULK_Si_PRIMITIVE
from unit.fixtures.nanoparticle import SI_NANOPARTICLE_SPHERE
from unit.utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "material_config, condition, expected_material_config",
    [
        (
            BULK_Si_PRIMITIVE,
            SphereCoordinateCondition(radius=6.0),
            SI_NANOPARTICLE_SPHERE,
        ),
    ],
)
def test_create_nanoparticle_by_condition(material_config, condition, expected_material_config):
    material = Material.create(material_config)
    nanoparticle = create_nanoparticle_from_material(
        material, condition=condition, center_around_atom=True, use_cartesian_coordinates=True
    )

    assert_two_entities_deep_almost_equal(nanoparticle, expected_material_config)


@pytest.mark.parametrize(
    "material_config, shape, parameters",
    [
        (BULK_Si_PRIMITIVE, NanoparticleShapesEnum.ICOSAHEDRON, {"noshells": 3}),
    ],
)
def test_create_nanoparticle_by_shape(material_config, shape, parameters):
    material = Material.create(material_config)
    nanoparticle = create_nanoparticle_by_shape(material, shape, parameters)
    assert isinstance(nanoparticle, Material)
