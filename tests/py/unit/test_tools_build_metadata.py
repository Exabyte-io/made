import pytest
from mat3ra.code.entity import InMemoryEntity, InMemoryEntityPydantic
from mat3ra.made.tools.build.metadata import BuildMetadata, MaterialMetadata
from pydantic import BaseModel


class ConfigWithToDict(InMemoryEntityPydantic):
    value: str = "test_dict"


# TODO: Remove this class when all configurations moved to Pydantic
class ConfigWithToJsonReturnsDict(BaseModel, InMemoryEntity):
    value: str = "test_json_dict"

    @property
    def _json(self):
        return {"value": self.value}


class ConfigWithToJsonReturnsStr(InMemoryEntityPydantic):
    value: str = "test_json_str"


class MockParameters(InMemoryEntityPydantic):
    param: str = "param_value"


def test_metadata_initialization_with_data():
    initial_data = {"existing_key": "existing_value", "build": [{"configuration": {"initial": "config"}}]}
    metadata = MaterialMetadata(**initial_data)

    assert metadata.existing_key == "existing_value"
    assert metadata.build[-1].configuration["initial"] == "config"


def test_metadata_empty_initialization():
    metadata = MaterialMetadata()
    assert metadata.model_dump(exclude_none=True) == {"build": [{"configuration": {}, "build_parameters": {}}]}


@pytest.mark.parametrize(
    "config_object, expected_dict",
    [
        (ConfigWithToDict(), {"value": "test_dict"}),
        (ConfigWithToJsonReturnsDict(), {"value": "test_json_dict"}),
        (ConfigWithToJsonReturnsStr(), {"value": "test_json_str"}),
    ],
)
def test_build_metadata_update(config_object, expected_dict):
    build_metadata = BuildMetadata()

    build_metadata.update(configuration=config_object)
    assert build_metadata.configuration == expected_dict

    build_metadata.configuration = {}
    build_metadata.update(build_parameters=config_object)
    assert build_metadata.build_parameters == expected_dict


def test_full_lifecycle_and_serialization():
    # Initialize with some existing data
    metadata = MaterialMetadata(existing_key="previous_material_metadata_value")

    # Update with a 'to_dict' config
    config = ConfigWithToDict(value="config_value")
    metadata.build[-1].update(configuration=config)

    params = MockParameters(param="param_value")
    metadata.build[-1].update(build_parameters=params)

    # Check the final state before serialization
    assert metadata.existing_key == "previous_material_metadata_value"
    assert metadata.build[-1].configuration["value"] == "config_value"
    assert metadata.build[-1].build_parameters["param"] == "param_value"

    # Check the final serialized dict
    final_dict = metadata.model_dump()
    expected_dict = {
        "existing_key": "previous_material_metadata_value",
        "build": [{"configuration": {"value": "config_value"}, "build_parameters": {"param": "param_value"}}],
    }
    assert final_dict == expected_dict
