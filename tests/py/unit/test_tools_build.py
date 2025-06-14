from mat3ra.made.material import Material
from mat3ra.made.tools.operations.core.binary import merge_materials
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.cuts import (
    CAVITY_MATERIAL_BASIS,
    FULL_MATERIAL,
    MERGED_CAVITY_SECTION_BASIS,
    MERGED_SECTION_CAVITY_BASIS,
    SECTION_MATERIAL_BASIS,
    SECTION_MATERIAL_BASIS_EXTRA_ATOM,
)
from mat3ra.made.tools.build import BaseConfigurationPydantic

section = Material.create({**FULL_MATERIAL, **SECTION_MATERIAL_BASIS})
cavity = Material.create({**FULL_MATERIAL, **CAVITY_MATERIAL_BASIS})
section_with_extra_atom = Material.create({**FULL_MATERIAL, **SECTION_MATERIAL_BASIS_EXTRA_ATOM})


def test_merge_materials():
    merged_material = merge_materials([section, cavity])
    merged_material_reverse = merge_materials([cavity, section])
    assertion_utils.assert_deep_almost_equal(merged_material.basis, MERGED_CAVITY_SECTION_BASIS)
    assertion_utils.assert_deep_almost_equal(merged_material_reverse.basis, MERGED_SECTION_CAVITY_BASIS)


def test_resolve_close_coordinates_basis():
    merged_material = merge_materials([section, cavity])
    merged_material_reverse = merge_materials([cavity, section])
    assertion_utils.assert_deep_almost_equal(merged_material.basis, MERGED_CAVITY_SECTION_BASIS)
    assertion_utils.assert_deep_almost_equal(merged_material_reverse.basis, MERGED_SECTION_CAVITY_BASIS)


def test_resolve_close_coordinates_basis_extra_atom():
    merged_material = merge_materials([section_with_extra_atom, cavity])
    merged_material_reverse = merge_materials([cavity, section_with_extra_atom])
    assertion_utils.assert_deep_almost_equal(merged_material.basis, MERGED_CAVITY_SECTION_BASIS)
    assertion_utils.assert_deep_almost_equal(merged_material_reverse.basis, MERGED_SECTION_CAVITY_BASIS)


def test_configuration():
    class MockConfiguration(BaseConfigurationPydantic):
        value_1: str = "test"
        value_2: int = 0

    class CompositeMockConfiguration(BaseConfigurationPydantic):
        component: MockConfiguration = MockConfiguration(value_1="component_value")

    mock_config = MockConfiguration(value_1="component_value", value_2=42)
    configuration = CompositeMockConfiguration(component=mock_config)
    assert isinstance(configuration.component, MockConfiguration)
    assert configuration.component.value_1 == "component_value"
    assert configuration.component.value_2 == 42

    config_dict = configuration.to_dict()
    assert config_dict["component"]["value_1"] == "component_value"
    assert config_dict["component"]["value_2"] == 42

    configuration_from_dict = CompositeMockConfiguration.from_dict(config_dict)
    assert isinstance(configuration_from_dict.component, MockConfiguration)
    assert configuration_from_dict.component.value_1 == "component_value"
    assert configuration_from_dict.component.value_2 == 42
