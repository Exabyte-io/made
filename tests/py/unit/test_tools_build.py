import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseConfigurationPydantic
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

section = Material.create({**FULL_MATERIAL, **SECTION_MATERIAL_BASIS})
cavity = Material.create({**FULL_MATERIAL, **CAVITY_MATERIAL_BASIS})
section_with_extra_atom = Material.create({**FULL_MATERIAL, **SECTION_MATERIAL_BASIS_EXTRA_ATOM})


@pytest.mark.parametrize(
    "materials_to_merge, expected_basis, expected_basis_reverse",
    [
        (
            [section, cavity],
            MERGED_CAVITY_SECTION_BASIS,
            MERGED_SECTION_CAVITY_BASIS,
        ),
        (
            [section_with_extra_atom, cavity],
            MERGED_CAVITY_SECTION_BASIS,
            MERGED_SECTION_CAVITY_BASIS,
        ),
    ],
)
def test_merge_materials(materials_to_merge, expected_basis, expected_basis_reverse):
    merged_material = merge_materials(materials_to_merge)
    merged_material_reverse = merge_materials(materials_to_merge[::-1])
    assertion_utils.assert_deep_almost_equal(merged_material.basis, expected_basis)
    assertion_utils.assert_deep_almost_equal(merged_material_reverse.basis, expected_basis_reverse)


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

    configuration_from_dict = CompositeMockConfiguration(**config_dict)
    assert isinstance(configuration_from_dict.component, MockConfiguration)
    assert configuration_from_dict.component.value_1 == "component_value"
    assert configuration_from_dict.component.value_2 == 42
