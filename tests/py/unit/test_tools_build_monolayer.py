import pytest
from mat3ra.made.material import Material

from mat3ra.made.tools.build.monolayer import create_monolayer
from .fixtures.bulk import BULK_Si_CONVENTIONAL, BULK_Si_PRIMITIVE
from .fixtures.generated.fixtures import BULK_GRAPHITE
from .utils import assert_two_entities_deep_almost_equal




@pytest.mark.parametrize(
    "material_config, vacuum, expected_vacuum_check",
    [
        (BULK_Si_PRIMITIVE, 10.0, True),
        (BULK_GRAPHITE, 15.0, True),
    ],
)
def test_create_monolayer(material_config, vacuum, expected_vacuum_check):
    """Test creating monolayer from different crystal types"""
    crystal = Material.create(material_config)
    monolayer = create_monolayer(crystal, vacuum=vacuum)
    
    # Basic checks
    assert isinstance(monolayer, Material)
    assert "Monolayer" in monolayer.name
    
    # Check vacuum was added correctly
    if expected_vacuum_check and vacuum > 0:
        assert monolayer.lattice.c > crystal.lattice.c
    
    # Check that we have fewer or equal atoms than the original (due to filtering)
    assert len(monolayer.basis.elements.values) <= len(crystal.basis.elements.values)


@pytest.mark.parametrize(
    "material_config, expected_lattice_type, expected_method, expected_miller_str",
    [
        (BULK_Si_CONVENTIONAL, "TRI", "_create_fcc_cub_monolayer", "111"),  # Si_CONVENTIONAL shows as TRI
        (BULK_GRAPHITE, "HEX", "_create_hex_monolayer", "001"),
    ],
)
def test_monolayer_builder_logic(material_config, expected_lattice_type, expected_method, expected_miller_str):
    """Test that the builder uses the correct method and Miller indices based on lattice type"""
    crystal = Material.create(material_config)
    
    # Convert enum to string for comparison
    lattice_type_str = crystal.lattice.type.value if hasattr(crystal.lattice.type, 'value') else str(crystal.lattice.type)
    assert lattice_type_str == expected_lattice_type
    
    configuration = MonolayerConfiguration.from_parameters(
        crystal=crystal,
        vacuum=10.0,
    )
    
    builder = MonolayerBuilder()
    monolayer = builder.get_material(configuration)
    
    # Check that we get a valid monolayer
    assert isinstance(monolayer, Material)
    assert "Monolayer" in monolayer.name
    assert f"({expected_miller_str})" in monolayer.name


def test_monolayer_configuration():
    """Test MonolayerConfiguration creation and serialization"""
    crystal = Material.create(BULK_Si_CONVENTIONAL)
    
    # Test from_parameters
    config = MonolayerConfiguration.from_parameters(
        crystal=crystal,
        vacuum=15.0,
    )
    
    assert config.crystal == crystal
    assert config.vacuum == 15.0
    
    # Test JSON serialization
    json_data = config._json
    assert "crystal" in json_data
    assert json_data["vacuum"] == 15.0


def test_monolayer_builder_name_update():
    """Test that material names are updated correctly with automatic Miller indices"""
    # Test FCC/CUB crystal (should get 111)
    crystal_si = Material.create(BULK_Si_CONVENTIONAL)
    original_name_si = crystal_si.name
    
    monolayer_si = create_monolayer(crystal_si, vacuum=10.0)
    
    assert monolayer_si.name == f"{original_name_si} - Monolayer (111)"
    assert monolayer_si.name != original_name_si
    
    # Test HEX crystal (should get 001)
    crystal_graphene = Material.create(GRAPHENE)
    original_name_graphene = crystal_graphene.name
    
    monolayer_graphene = create_monolayer(crystal_graphene, vacuum=10.0)
    
    assert monolayer_graphene.name == f"{original_name_graphene} - Monolayer (001)"
    assert monolayer_graphene.name != original_name_graphene


def test_automatic_miller_indices():
    """Test that Miller indices are automatically determined correctly"""
    # HEX crystal should use (0,0,1)
    graphene = Material.create(GRAPHENE)
    monolayer_graphene = create_monolayer(graphene, vacuum=10.0)
    assert "(001)" in monolayer_graphene.name
    
    # FCC/CUB crystal should use (1,1,1)  
    silicon = Material.create(BULK_Si_CONVENTIONAL)
    monolayer_silicon = create_monolayer(silicon, vacuum=10.0)
    assert "(111)" in monolayer_silicon.name 
