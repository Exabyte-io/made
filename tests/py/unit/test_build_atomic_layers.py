from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer

from mat3ra.made.tools.build.slab.helpers import select_slab_termination, get_slab_terminations
from mat3ra.standata.materials import Materials

from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.builders import AtomicLayersUniqueRepeatedBuilder
from mat3ra.made.tools.build.slab.configuration import (
    CrystalLatticePlanesConfiguration,
    AtomicLayersUniqueRepeatedConfiguration,
)
from mat3ra.made.tools.modify import wrap_to_unit_cell
from mat3ra.made.tools.operations.core.unary import translate
from mat3ra.made.tools.analyze.other import get_surface_atom_indices, get_local_extremum_atom_index
from mat3ra.made.tools.enums import SurfaceTypes
from unit.utils import assert_two_entities_deep_almost_equal


SrTiO_BULK = Material.create(Materials.get_by_name_first_match("SrTiO3"))
material = SrTiO_BULK
MILLER_INDICES = (1, 1, 0)
NUMBER_OF_LAYERS = 3


def process_termination(material, termination):
    atomic_layers_repeated_config = AtomicLayersUniqueRepeatedConfiguration(
        crystal=material,
        miller_indices=MILLER_INDICES,
        termination_top=termination,
        number_of_repetitions=NUMBER_OF_LAYERS,
    )

    atomic_layers_repeated_orthogonal_c = AtomicLayersUniqueRepeatedBuilder().get_material(
        atomic_layers_repeated_config
    )
    translation_vector = atomic_layers_repeated_config.get_translation_vector(termination)
    atomic_layers_repeated_terminated = translate(atomic_layers_repeated_orthogonal_c, translation_vector)
    slab = wrap_to_unit_cell(atomic_layers_repeated_terminated)
    return slab


def get_topmost_atom_element(slab):
    center_coord = [0.5, 0.5, 0.5]
    topmost_index = get_local_extremum_atom_index(slab, center_coord, "max", vicinity=1000.0)

    return slab.basis.elements.get_element_value_by_index(topmost_index)


def test_termination_translation():
    material = SrTiO_BULK

    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
        material=material, miller_indices=MILLER_INDICES
    )
    terminations = crystal_lattice_planes_analyzer.terminations

    termination = select_slab_termination(terminations, "SrTiO")
    slab_1 = process_termination(material, termination)
    topmost_atom_element = get_topmost_atom_element(slab_1)
    assert topmost_atom_element == "Sr"

    termination = select_slab_termination(terminations, "O2")
    slab_2 = process_termination(material, termination)

    topmost_atom_element = get_topmost_atom_element(slab_2)
    assert topmost_atom_element == "O"
