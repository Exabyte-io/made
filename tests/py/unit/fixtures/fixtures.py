from ase.build import bulk
from mat3ra.made.material import Material
from mat3ra.made.tools.build.interface.termination_pair import TerminationPair
from mat3ra.made.tools.build.slab import SlabConfiguration, create_slab, get_terminations
from mat3ra.made.tools.convert import from_ase
from pymatgen.analysis.elasticity.strain import Strain
from pymatgen.core.interface import Interface
from unit.utils import atoms_to_interface_structure

# ASE Atoms fixtures
substrate = bulk("Si", cubic=True)
film = bulk("Cu", cubic=True)

# TODO: move to interface.py
INTERFACE_ATOMS = substrate + film
INTERFACE_ATOMS.set_tags([1] * len(substrate) + [2] * len(film))

# Material fixtures
SUBSTRATE_MATERIAL = Material.create(from_ase(substrate))
FILM_MATERIAL = Material.create(from_ase(film))

SUBSTRATE_CONFIGURATION = SlabConfiguration(bulk=SUBSTRATE_MATERIAL, thickness=3)
FILM_CONFIGURATION = SlabConfiguration(bulk=FILM_MATERIAL)

substrate_terminations = get_terminations(SUBSTRATE_CONFIGURATION)
film_terminations = get_terminations(FILM_CONFIGURATION)

# Pymatgen Interface fixtures
INTERFACE_TERMINATION_PAIR: TerminationPair = TerminationPair(
    film_terminations[0],
    substrate_terminations[0],
)
INTERFACE_TERMINATION_AS_STR = str(INTERFACE_TERMINATION_PAIR)

interface_structure = atoms_to_interface_structure(INTERFACE_ATOMS)
interface_dict = interface_structure.as_dict()

INTERFACE_STRUCTURE = Interface.from_dict(interface_dict)
# Create properties that are assigned during interface creation in ZSL algorithm as dict and json
INTERFACE_PROPERTIES_MOCK = {
    "film_transformation": [[2.0, 0.0], [0.0, 2.0]],
    "substrate_transformation": [[1.0, 0.0], [0.0, 1.0]],
    "strain": Strain([[0.004746364, -0.0, -0.0], [-0.0, 0.004746364, 0.0], [-0.0, 0.0, -0.0]]),
    "von_mises_strain": 0.001,
    "mean_abs_strain": 0.00105,
    "termination": (INTERFACE_TERMINATION_PAIR.film_termination, INTERFACE_TERMINATION_PAIR.substrate_termination),
}
INTERFACE_PROPERTIES_JSON = {
    "film_transformation": [[2.0, 0.0], [0.0, 2.0]],
    "substrate_transformation": [[1.0, 0.0], [0.0, 1.0]],
    "strain": [[0.004746364, -0.0, -0.0], [-0.0, 0.004746364, 0.0], [-0.0, 0.0, -0.0]],
    "von_mises_strain": 0.001,
    "termination": INTERFACE_TERMINATION_AS_STR,
    "mean_abs_strain": 0.00105,
}

# Add properties to interface structure
INTERFACE_STRUCTURE.interface_properties = INTERFACE_PROPERTIES_MOCK
INTERFACE_NAME = "Cu4(001)-Si8(001), Interface, Strain 0.062pct"


clean_material = Material.create_default()
slab_111_config = SlabConfiguration(
    bulk=clean_material,
    miller_indices=(1, 1, 1),
    thickness=4,
    vacuum=6,
    xy_supercell_matrix=[[1, 0], [0, 1]],
    use_orthogonal_z=True,
)
t_111 = get_terminations(slab_111_config)[0]
SLAB_111 = create_slab(slab_111_config, t_111)

slab_001_config = SlabConfiguration(
    bulk=clean_material,
    miller_indices=(0, 0, 1),
    thickness=3,
    vacuum=3,
    xy_supercell_matrix=[[2, 0], [0, 1]],
    use_orthogonal_z=True,
)
t_001 = get_terminations(slab_001_config)[0]
SLAB_001 = create_slab(slab_001_config, t_001)
