from mat3ra.made.tools.build.interface import InterfaceDataHolder, patch_interface_with_mean_abs_strain
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.interface import Interface


# Create a mock interface structure
nickel_lattice = Lattice.hexagonal(a=2.49, c=4.00)
nickel = Structure(nickel_lattice, ["Ni", "Ni"], [[1 / 3, 2 / 3, 0.5], [2 / 3, 1 / 3, 0.5]])

slab_gen = SlabGenerator(nickel, miller_index=(1, 1, 1), min_slab_size=10, min_vacuum_size=10)
nickel_slab = slab_gen.get_slab()

graphene_lattice = Lattice.hexagonal(a=2.46, c=10)
graphene = Structure(graphene_lattice, ["C", "C"], [[0, 0, 0.1], [1 / 3, 2 / 3, 0.1]])

slab_gen = SlabGenerator(graphene, miller_index=(0, 0, 1), min_slab_size=10, min_vacuum_size=10)
graphene_slab = slab_gen.get_slab()

interface_structure = Interface.from_slabs(
    substrate_slab=nickel_slab,
    film_slab=graphene_slab,
    gap=2.0,
    vacuum_over_film=10.0,
    interface_properties={"termination": "top", "orientation": "aligned"},
)

interface_structure.interface_properties["strain"] = 0.1
patch_interface_with_mean_abs_strain(interface_structure)


def test_add_data_entries():
    interfaces_data = InterfaceDataHolder()
    interfaces_data.add_data_entries([interface_structure])
    assert len(interfaces_data.get_interfaces_for_termination(0)) == 1


def test_get_interfaces_for_termination():
    interfaces_data = InterfaceDataHolder()
    interfaces_data.add_data_entries([interface_structure])
    assert interfaces_data.get_interfaces_for_termination(0)[0] == interface_structure
