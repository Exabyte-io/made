from ase.build import bulk
from mat3ra.made.material import Material
from mat3ra.made.tools.build.interface import get_optimal_film_displacement
from mat3ra.made.tools.convert import from_ase
from mat3ra.made.tools.convert.utils import InterfacePartsEnum
from mat3ra.made.tools.modify import (
    add_vacuum,
    filter_by_circle_projection,
    filter_by_ids,
    filter_by_label,
    filter_by_layers,
    filter_by_rectangle_projection,
    filter_by_sphere,
    filter_by_triangle_projection,
    interface_displace_part,
    remove_vacuum,
    rotate,
    translate_to_z_level,
)
from mat3ra.utils import assertion as assertion_utils

from .fixtures.cell import SI_CONVENTIONAL_CELL, SI_CONVENTIONAL_CELL_FILTERED
from .fixtures.interface import GRAPHENE_NICKEL_INTERFACE
from .fixtures.slab import SI_SLAB_001, SI_SLAB_001_WITH_VACUUM
from .utils import assert_two_entities_deep_almost_equal

COMMON_PART = {
    "units": "crystal",
    "labels": [],
}

expected_basis_layers_section = {
    "elements": [
        {"id": 0, "value": "Si"},
        {"id": 3, "value": "Si"},
        {"id": 5, "value": "Si"},
        {"id": 6, "value": "Si"},
    ],
    "coordinates": [
        {"id": 0, "value": [0.5, 0.0, 0.0]},
        {"id": 3, "value": [0.25, 0.75, 0.25]},
        {"id": 5, "value": [0.75, 0.25, 0.25]},
        {"id": 6, "value": [0.0, 0.5, 0.0]},
    ],
    **COMMON_PART,
}

expected_basis_layers_cavity = {
    "elements": [
        {"id": 1, "value": "Si"},
        {"id": 2, "value": "Si"},
        {"id": 4, "value": "Si"},
        {"id": 7, "value": "Si"},
    ],
    "coordinates": [
        {"id": 1, "value": [0.25, 0.25, 0.75]},
        {"id": 2, "value": [0.5, 0.5, 0.5]},
        {"id": 4, "value": [0.0, 0.0, 0.5]},
        {"id": 7, "value": [0.75, 0.75, 0.75]},
    ],
    **COMMON_PART,
}


expected_basis_sphere_cluster = {
    "elements": [{"id": 2, "value": "Si"}],
    "coordinates": [{"id": 2, "value": [0.5, 0.5, 0.5]}],
    **COMMON_PART,
}

expected_basis_sphere_cavity = {
    "elements": [
        {"id": 0, "value": "Si"},
        {"id": 1, "value": "Si"},
        {"id": 3, "value": "Si"},
        {"id": 4, "value": "Si"},
        {"id": 5, "value": "Si"},
        {"id": 6, "value": "Si"},
        {"id": 7, "value": "Si"},
    ],
    "coordinates": [
        {"id": 0, "value": [0.5, 0.0, 0.0]},
        {"id": 1, "value": [0.25, 0.25, 0.75]},
        {"id": 3, "value": [0.25, 0.75, 0.25]},
        {"id": 4, "value": [0.0, 0.0, 0.5]},
        {"id": 5, "value": [0.75, 0.25, 0.25]},
        {"id": 6, "value": [0.0, 0.5, 0.0]},
        {"id": 7, "value": [0.75, 0.75, 0.75]},
    ],
    **COMMON_PART,
}

CRYSTAL_RADIUS = 0.25  # in crystal coordinates
CRYSTAL_CENTER_3D = [0.5, 0.5, 0.5]  # in crystal coordinates


def test_filter_by_ids():
    material = Material.create(SI_CONVENTIONAL_CELL)
    material_filtered = filter_by_ids(material, [0, 2])
    # Move to fixtures
    expected_material = Material.create(SI_CONVENTIONAL_CELL_FILTERED)
    assert_two_entities_deep_almost_equal(material_filtered.basis, expected_material.basis)


def test_filter_by_label():
    # TODO: create in fixtures
    substrate = bulk("Si", cubic=True)
    film = bulk("Cu", cubic=True)
    interface_atoms = substrate + film
    interface_atoms.set_tags([1] * len(substrate) + [2] * len(film))
    material_interface = Material.create(from_ase(interface_atoms))
    film_extracted = filter_by_label(material_interface, 2)
    film_material = Material.create(from_ase(film))

    # Ids of filtered elements will be missing, comparing the resulting values
    assert [el for el in film_material.basis.elements.values] == [el for el in film_extracted.basis.elements.values]


def test_filter_by_layers():
    material = Material.create(SI_CONVENTIONAL_CELL)
    section = filter_by_layers(material=material, central_atom_id=0, layer_thickness=3.0)
    cavity = filter_by_layers(material=material, central_atom_id=0, layer_thickness=3.0, invert_selection=True)
    assert_two_entities_deep_almost_equal(section.basis, expected_basis_layers_section)
    assert_two_entities_deep_almost_equal(cavity.basis, expected_basis_layers_cavity)


def test_filter_by_sphere():
    material = Material.create(SI_CONVENTIONAL_CELL)
    cluster = filter_by_sphere(material, center_coordinate=CRYSTAL_CENTER_3D, radius=CRYSTAL_RADIUS)
    cavity = filter_by_sphere(material, center_coordinate=CRYSTAL_CENTER_3D, radius=CRYSTAL_RADIUS, invert=True)
    assert_two_entities_deep_almost_equal(cluster.basis, expected_basis_sphere_cluster)
    assert_two_entities_deep_almost_equal(cavity.basis, expected_basis_sphere_cavity)


def test_filter_by_circle_projection():
    material = Material.create(SI_CONVENTIONAL_CELL)
    # Small cylinder in the middle of the cell containing the central atom will be removed -- the same as with sphere
    section = filter_by_circle_projection(material, 0.5, 0.5, CRYSTAL_RADIUS)
    cavity = filter_by_circle_projection(material, 0.5, 0.5, CRYSTAL_RADIUS, invert_selection=True)
    assert_two_entities_deep_almost_equal(section.basis, expected_basis_sphere_cluster)
    assert_two_entities_deep_almost_equal(cavity.basis, expected_basis_sphere_cavity)


def test_filter_by_rectangle_projection():
    material = Material.create(SI_CONVENTIONAL_CELL)
    # Default will contain all the atoms
    section = filter_by_rectangle_projection(material)
    assert_two_entities_deep_almost_equal(material.basis, section.basis)


def test_filter_by_triangle_projection():
    # Small prism in the middle of the cell containing the central atom will be removed -- the same as with sphere
    material = Material.create(SI_CONVENTIONAL_CELL)
    section = filter_by_triangle_projection(material, [0.4, 0.4], [0.4, 0.5], [0.5, 0.5])
    cavity = filter_by_triangle_projection(material, [0.4, 0.4], [0.4, 0.5], [0.5, 0.5], invert_selection=True)
    assert_two_entities_deep_almost_equal(section.basis, expected_basis_sphere_cluster)
    assert_two_entities_deep_almost_equal(cavity.basis, expected_basis_sphere_cavity)


def test_add_vacuum():
    material = Material.create(SI_SLAB_001)
    material_with_vacuum = add_vacuum(material, 5.0)
    assert_two_entities_deep_almost_equal(material_with_vacuum, SI_SLAB_001_WITH_VACUUM)


def test_remove_vacuum():
    material_with_vacuum = Material.create(SI_SLAB_001_WITH_VACUUM)
    vacuum = 7.368
    material_with_no_vacuum = remove_vacuum(material_with_vacuum, from_top=True, from_bottom=True, fixed_padding=0)
    expected_material = add_vacuum(material_with_no_vacuum, vacuum)
    # to compare correctly, we need to translate the expected material to the bottom
    # as it is effectively moved down when vacuum is removed (set to 0), use atol=1e-3 to account for the translation
    reference_material = Material.create(SI_SLAB_001)
    reference_material_translated_down = translate_to_z_level(reference_material, z_level="bottom")
    # compare absolute coordinates

    assert_two_entities_deep_almost_equal(expected_material, reference_material_translated_down, atol=1e-3)


def test_rotate():
    material = Material.create(SI_SLAB_001)
    # Rotation around Z and X axis will be equivalent for the original material for hist basis in terms of coordinates
    rotated_material = rotate(material, [0, 0, 1], 180)
    rotated_material = rotate(rotated_material, [1, 0, 0], 180)
    assertion_utils.assert_deep_almost_equal(
        material.basis.coordinates.values.sort(), rotated_material.basis.coordinates.values.sort()
    )
    assertion_utils.assert_deep_almost_equal(material.lattice, rotated_material.lattice)


def test_displace_interface():
    material = Material.create(GRAPHENE_NICKEL_INTERFACE)
    expected_coordinates = [
        {"id": 0, "value": [0.666666667, 0.666666667, 0.350869517]},
        {"id": 1, "value": [-0.0, 0.0, 0.425701769]},
        {"id": 2, "value": [0.333333333, 0.333333333, 0.500534021]},
        {"id": 3, "value": [0.433333333, 0.533333333, 0.911447347]},
        {"id": 4, "value": [0.766666667, 0.866666667, 0.911447347]},
    ]
    expected_labels = GRAPHENE_NICKEL_INTERFACE["basis"]["labels"]
    displaced_material = interface_displace_part(
        material, [0.1, 0.2, 0.3], InterfacePartsEnum.FILM, use_cartesian_coordinates=False
    )
    assertion_utils.assert_deep_almost_equal(expected_coordinates, displaced_material.basis.coordinates.to_dict())
    assertion_utils.assert_deep_almost_equal(expected_labels, displaced_material.basis.labels.to_dict())


def test_displace_interface_optimized():
    material = Material.create(GRAPHENE_NICKEL_INTERFACE)
    expected_coordinates = [
        {"id": 0, "value": [0.666666667, 0.666666667, 0.350869517]},
        {"id": 1, "value": [-0.0, 0.0, 0.425701769]},
        {"id": 2, "value": [0.333333333, 0.333333333, 0.500534021]},
        {"id": 3, "value": [0.285973954, 0.203945038, 0.611447347]},
        {"id": 4, "value": [0.619307288, 0.537278372, 0.611447347]},
    ]
    expected_labels = GRAPHENE_NICKEL_INTERFACE["basis"]["labels"]

    optimal_displacement = get_optimal_film_displacement(
        material, grid_size_xy=(10, 10), grid_range_x=(-0.5, 0.5), grid_range_y=(-0.5, 0.5)
    )
    displaced_material = interface_displace_part(material, optimal_displacement, use_cartesian_coordinates=True)
    assertion_utils.assert_deep_almost_equal(expected_coordinates, displaced_material.basis.coordinates.to_dict())
    assertion_utils.assert_deep_almost_equal(expected_labels, displaced_material.basis.labels.to_dict())
