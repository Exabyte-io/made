from ase.build import bulk
from mat3ra.made.material import Material
from mat3ra.made.tools.build.utils import merge_materials
from mat3ra.made.tools.convert import from_ase
from mat3ra.made.tools.modify import filter_by_layers
from mat3ra.utils import assertion as assertion_utils

ase_ni = bulk("Ni", "fcc", a=3.52, cubic=True)
material = Material(from_ase(ase_ni))
section = filter_by_layers(material, central_atom_id=0, layer_thickness=1.0)
cavity = filter_by_layers(material, central_atom_id=0, layer_thickness=1.0, invert_selection=True)

# Change 0th element
section.basis.elements.values[0] = "Ge"

# Add element to cavity for collision test
cavity.basis.elements.add_item("S", id=4)
coordinate_value = section.basis.coordinates.values[1]
cavity.basis.coordinates.add_item(coordinate_value, id=4)

expected_merged_material_basis = {
    "elements": [{"id": 0, "value": "Ge"}, {"id": 1, "value": "Ni"}, {"id": 2, "value": "Ni"}, {"id": 4, "value": "S"}],
    "coordinates": [
        {"id": 0, "value": [0.0, 0.0, 0.0]},
        {"id": 1, "value": [0.0, 0.5, 0.5]},
        {"id": 2, "value": [0.5, 0.0, 0.5]},
        {"id": 4, "value": [0.5, 0.5, 0.0]},
    ],
    "labels": [],
}


expected_merged_material_reverse_basis = {
    "elements": [
        {"id": 1, "value": "Ni"},
        {"id": 2, "value": "Ni"},
        {"id": 0, "value": "Ge"},
        {"id": 3, "value": "Ni"},
    ],
    "coordinates": [
        {"id": 1, "value": [0.0, 0.5, 0.5]},
        {"id": 2, "value": [0.5, 0.0, 0.5]},
        {"id": 0, "value": [0.0, 0.0, 0.0]},
        {"id": 3, "value": [0.5, 0.5, 0.0]},
    ],
    "labels": [],
}


def test_merge_materials():
    merged_material = merge_materials([section, cavity])
    merged_material_reverse = merge_materials([cavity, section])
    assertion_utils.assert_deep_almost_equal(merged_material.basis, expected_merged_material_basis)
    assertion_utils.assert_deep_almost_equal(merged_material_reverse.basis, expected_merged_material_reverse_basis)
