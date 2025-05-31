from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.build.utils import stack_two_materials


def test_stack_two_materials_z_direction():
    """Test stacking two materials along z direction"""
    material1 = Material.create_default()
    material2 = Material.create_default()

    original_c = material1.lattice.c
    original_atom_count = len(material1.basis.coordinates.values)

    stacked_material_z = stack_two_materials(material1, material2, AxisEnum.z)

    # Check that c lattice parameter is doubled
    assert abs(stacked_material_z.lattice.c - 2 * original_c) < 1e-10
    # Check that a and b parameters remain the same
    assert abs(stacked_material_z.lattice.a - material1.lattice.a) < 1e-10
    assert abs(stacked_material_z.lattice.b - material1.lattice.b) < 1e-10
    # Check that atom count is doubled
    assert len(stacked_material_z.basis.coordinates.values) == 2 * original_atom_count


def test_stack_two_materials_x_direction():
    """Test stacking two materials along x direction"""
    material1 = Material.create_default()
    material2 = Material.create_default()

    original_a = material1.lattice.a
    original_atom_count = len(material1.basis.coordinates.values)

    stacked_material_x = stack_two_materials(material1, material2, AxisEnum.x)

    # Check that a lattice parameter is doubled
    assert abs(stacked_material_x.lattice.a - 2 * original_a) < 1e-10
    # Check that b and c parameters remain the same
    assert abs(stacked_material_x.lattice.b - material1.lattice.b) < 1e-10
    assert abs(stacked_material_x.lattice.c - material1.lattice.c) < 1e-10
    # Check that atom count is doubled
    assert len(stacked_material_x.basis.coordinates.values) == 2 * original_atom_count
