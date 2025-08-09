from types import SimpleNamespace

import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.helpers import StackComponentDict, create_heterostructure

from .fixtures.bulk import BULK_Hf2O_MCL, BULK_Si_CONVENTIONAL, BULK_SiO2, BULK_TiN

PRECISION = 1e-3

Si_SiO2_Hf2O_HETEROSTRUCTURE_TEST_CASE = (
    [
        SimpleNamespace(bulk_config=BULK_Si_CONVENTIONAL, miller_indices=(0, 0, 1), number_of_layers=4),
        SimpleNamespace(bulk_config=BULK_SiO2, miller_indices=(1, 1, 1), number_of_layers=5),
        SimpleNamespace(bulk_config=BULK_Hf2O_MCL, miller_indices=(0, 0, 1), number_of_layers=2),
        SimpleNamespace(bulk_config=BULK_TiN, miller_indices=(1, 1, 1), number_of_layers=5),
    ],
    [1.5, 1.0, 1.0],  # gaps
    10.0,  # vacuum
)


@pytest.mark.parametrize("layers, gaps, vacuum", [Si_SiO2_Hf2O_HETEROSTRUCTURE_TEST_CASE])
def test_create_heterostructure_simple(layers, gaps, vacuum):
    # Convert raw test data to Pydantic models
    stack_components = []
    for layer in layers:
        component_data = {
            "crystal": Material.create(layer.bulk_config),
            "miller_indices": layer.miller_indices,
            "thickness": layer.number_of_layers,
        }
        component = StackComponentDict(**component_data)
        stack_components.append(component)

    heterostructure = create_heterostructure(
        stack_component_dicts=stack_components,
        gaps=gaps,
        vacuum=vacuum,
    )

    assert isinstance(heterostructure, Material)
    elements = set()
    for layer in layers:
        elements.update(Material.create(layer.bulk_config).basis.elements.values)
    assert set(heterostructure.basis.elements.values) == elements

    assert len(heterostructure.basis.elements.values) > 0
