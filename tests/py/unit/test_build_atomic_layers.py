from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.analyze.other import get_local_extremum_atom_index
from mat3ra.made.tools.build.slab.builders import AtomicLayersUniqueRepeatedBuilder, CrystalLatticePlanesBuilder
from mat3ra.made.tools.build.slab.configuration import (
    AtomicLayersUniqueRepeatedConfiguration,
    CrystalLatticePlanesConfiguration,
)
from mat3ra.made.tools.build.slab.helpers import select_slab_termination
from unit.fixtures.generated.fixtures import SrTiO3_BULK_MATERIAL

NUMBER_OF_LAYERS = 1


def process_termination(material, miller_indices, termination):
    crystal_lattice_planes_config = CrystalLatticePlanesConfiguration(crystal=material, miller_indices=miller_indices)
    crystal_lattice_planes_material = CrystalLatticePlanesBuilder().get_material(crystal_lattice_planes_config)

    atomic_layers_config = AtomicLayersUniqueRepeatedConfiguration(
        crystal=crystal_lattice_planes_material,
        miller_indices=miller_indices,
        termination_top=termination,
        number_of_repetitions=NUMBER_OF_LAYERS,
    )

    atomic_layers_builder = AtomicLayersUniqueRepeatedBuilder()
    atomic_layers_material = atomic_layers_builder.get_material(atomic_layers_config)

    return atomic_layers_material


def get_topmost_atom_element(slab):
    center_coord = [0.5, 0.5, 0.5]
    topmost_index = get_local_extremum_atom_index(slab, center_coord, "max", vicinity=1000.0)
    return slab.basis.elements.get_element_value_by_index(topmost_index)


def test_termination_translation_srtio3():
    material = SrTiO3_BULK_MATERIAL
    miller_indices = (1, 1, 0)

    analyzer = CrystalLatticePlanesMaterialAnalyzer(material=material, miller_indices=miller_indices)
    terminations = analyzer.terminations

    termination_1 = select_slab_termination(terminations, "SrTiO")
    slab_1 = process_termination(material, miller_indices, termination_1)
    topmost_atom_element_1 = get_topmost_atom_element(slab_1)

    termination_2 = select_slab_termination(terminations, "O2")
    slab_2 = process_termination(material, miller_indices, termination_2)
    topmost_atom_element_2 = get_topmost_atom_element(slab_2)

    assert topmost_atom_element_1 == "Sr"
    assert topmost_atom_element_2 == "O"
