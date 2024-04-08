from mat3ra.made.material import Material

REFERENCE_OBJECT_1 = {"key1": "value1", "key2": "value2"}


def test_create():
    material = Material.create(Material.default_config)
    assert material.to_json() == Material.default_config
    assert material.name == Material.default_config["name"]
