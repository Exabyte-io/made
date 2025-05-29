from typing import List, Union

from mat3ra.esse.models.material.reusable.two_dimensional.atomic_layers import AtomicLayersSchema
from mat3ra.esse.models.material.reusable.two_dimensional.atomic_layers_unique import (
    AtomicLayersUniqueSchema,
)
from mat3ra.esse.models.material.reusable.two_dimensional.atomic_layers_unique_repeated import (
    AtomicLayersUniqueRepeatedSchema,
)
from mat3ra.esse.models.material.reusable.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import (
    SlabConfigurationSchema,
    AxisEnum,
    VacuumConfigurationSchema,
)
from pymatgen.core import Lattice

from mat3ra.made.material import Material
from .termination import Termination
from .. import BaseConfigurationPydantic
from ...convert import to_pymatgen, from_pymatgen
from ...third_party import PymatgenSlab
from ...third_party import PymatgenSlabGenerator, label_pymatgen_slab_termination
from ...third_party import PymatgenSpacegroupAnalyzer


class CrystalLatticePlanes(CrystalLatticePlanesSchema):

    crystal: Material
    _pymatgen_slabs: List[PymatgenSlab]

    def __init__(self, **kwargs):
        crystal = kwargs.pop("crystal", None) or Material.create_default()
        use_conventional = kwargs.get("use_conventional_cell", False)
        if use_conventional:
            # TODO: use from LatticeAnalyzer
            crystal = Material.create(
                from_pymatgen(PymatgenSpacegroupAnalyzer(to_pymatgen(crystal)).get_conventional_standard_structure())
            )
        kwargs["crystal"] = crystal
        super().__init__(**kwargs)
        self._pymatgen_slabs = self._generate_pymatgen_slabs()

    def _generate_pymatgen_slabs(
        self,
        min_slab_size=1,
        min_vacuum_size=1,
        in_unit_planes: bool = True,
        make_primitive: bool = False,
        symmetrize: bool = False,
    ) -> List[PymatgenSlab]:  # type: ignore
        generator = PymatgenSlabGenerator(
            initial_structure=to_pymatgen(self.crystal),
            miller_index=self.miller_indices,
            min_slab_size=min_slab_size,
            min_vacuum_size=min_vacuum_size,
            in_unit_planes=in_unit_planes,
            primitive=make_primitive,
        )
        raw_slabs = generator.get_slabs(
            # We need to preserve symmetric slabs for different terminations at the surface
            symmetrize=symmetrize
        )

        return raw_slabs

    def get_terminations(self):
        return [Termination.from_string(label_pymatgen_slab_termination(slab)) for slab in self._pymatgen_slabs]

    def get_slabs(
        self,
        min_slab_size=1,
        min_vacuum_size=1,
        in_unit_planes: bool = True,
        make_primitive: bool = False,
        symmetrize: bool = False,
    ) -> List[Material]:

        pymatgen_slabs = self._generate_pymatgen_slabs(
            min_slab_size=min_slab_size,
            min_vacuum_size=min_vacuum_size,
            in_unit_planes=in_unit_planes,
            make_primitive=make_primitive,
            symmetrize=symmetrize,
        )
        return [Material.create(from_pymatgen(slab)) for slab in pymatgen_slabs]

    @staticmethod
    def _strip_vacuum_from_pymatgen_slab(slab: PymatgenSlab, vacuum_to_remove: float) -> PymatgenSlab:
        """Return a copy of *slab* with its c lattice vector shortened by
        *vacuum_to_remove* Å (negative values add vacuum).  Works only when the
        c‑vector is the out‑of‑plane direction, as is the case for PymatgenSlab."""
        old_c = slab.lattice.c
        new_c = old_c - vacuum_to_remove
        scale = old_c / new_c

        # new lattice
        a, b, c_vec = slab.lattice.matrix
        new_lattice = Lattice([a, b, c_vec * (new_c / old_c)])

        # rescale fractional z so atoms stay in place
        new_frac = slab.frac_coords.copy()
        new_frac[:, 2] *= scale
        new_frac %= 1.0  # wrap back into cell

        return PymatgenSlab(
            lattice=new_lattice,
            species=slab.species_and_occu,
            coords=new_frac,
            miller_index=slab.miller_index,
            oriented_unit_cell=slab.oriented_unit_cell,
            shift=slab.shift,
            scale_factor=slab.scale_factor,
            reorient_lattice=slab.reorient_lattice,
            site_properties=slab.site_properties,
            energy=slab.energy,
        )


class AtomicLayers(AtomicLayersSchema, CrystalLatticePlanes):
    pass


class AtomicLayersUnique(AtomicLayersUniqueSchema, CrystalLatticePlanes):
    pass


class AtomicLayersUniqueRepeated(AtomicLayersUniqueRepeatedSchema, CrystalLatticePlanes):
    crystal: Material
    termination_top: Termination

    def _generate_pymatgen_slabs(self, **kwargs) -> List[PymatgenSlab]:
        raw_slabs = super()._generate_pymatgen_slabs(
            min_slab_size=self.number_of_repetitions,
            min_vacuum_size=kwargs.get("min_vacuum_size", 1),
            in_unit_planes=kwargs.get("in_unit_planes", True),
            make_primitive=kwargs.get("make_primitive", False),
            symmetrize=kwargs.get("symmetrize", False),
        )

        min_vacuum_size = kwargs.get("min_vacuum_size", 1)
        in_unit_planes = kwargs.get("in_unit_planes", True)
        added_vacuum = min_vacuum_size * self.crystal.lattice.c if in_unit_planes else min_vacuum_size

        no_vac_slabs = [self._strip_vacuum_from_pymatgen_slab(slab, added_vacuum) for slab in raw_slabs]
        return no_vac_slabs


class VacuumConfiguration(VacuumConfigurationSchema):
    """
    Configuration for adding vacuum to a material.

    Args:
        direction (AxisEnum): The direction to add vacuum along (x, y, or z).
        size (float): Size of the vacuum gap in angstroms. Defaults to 10.0.
    """

    direction: AxisEnum = VacuumConfigurationSchema.model_fields["direction"].default
    size: float = VacuumConfigurationSchema.model_fields["size"].default


class SlabConfiguration(SlabConfigurationSchema, BaseConfigurationPydantic):
    type: str = "SlabConfiguration"
    stack_components: List[Union[AtomicLayersUniqueRepeated, VacuumConfiguration]]
    xy_supercell_matrix: List[List[int]] = SlabConfigurationSchema.model_fields["xy_supercell_matrix"].default
    direction: AxisEnum = AxisEnum.z

    @property
    def bulk(self) -> Material:
        atomic_layers = self.stack_components[0]
        if hasattr(atomic_layers, "crystal"):
            return atomic_layers.crystal
        else:
            raise AttributeError("No crystal found in stack components")

    @property
    def miller_indices(self) -> tuple:
        atomic_layers = self.stack_components[0]
        if hasattr(atomic_layers, "miller_indices"):
            return atomic_layers.miller_indices
        else:
            raise AttributeError("No miller_indices found in stack components")

    @property
    def number_of_layers(self) -> int:
        atomic_layers = self.stack_components[0]
        if hasattr(atomic_layers, "number_of_repetitions"):
            return atomic_layers.number_of_repetitions
        else:
            raise AttributeError("No number_of_repetitions found in stack components")

    @classmethod
    def from_parameters(
        cls,
        bulk: Material,
        miller_indices: tuple = (0, 0, 1),
        number_of_layers: int = 1,
        vacuum: float = 10.0,
        xy_supercell_matrix: List[List[int]] = None,
        use_conventional_cell: bool = False,
        termination_index: int = 0,
        **kwargs,
    ) -> "SlabConfiguration":
        """
        Create SlabConfiguration from parameters.

        Args:
            bulk: The bulk material
            miller_indices: Miller indices for the surface
            number_of_layers: Number of atomic layers
            vacuum: Vacuum size in Angstroms
            xy_supercell_matrix: Supercell matrix for xy directions
            use_conventional_cell: Whether to use conventional cell
            termination_index: Index of the termination to use
        """
        if xy_supercell_matrix is None:
            xy_supercell_matrix = [[1, 0], [0, 1]]

        temp_planes = CrystalLatticePlanes(
            crystal=bulk,
            miller_indices=miller_indices,
            use_conventional_cell=use_conventional_cell,
        )

        available_terminations = temp_planes.get_terminations()
        termination = available_terminations[termination_index]

        atomic_layers = AtomicLayersUniqueRepeated(
            crystal=bulk,
            miller_indices=miller_indices,
            number_of_repetitions=number_of_layers,
            use_conventional_cell=use_conventional_cell,
            termination_top=termination,
        )

        vacuum_config = VacuumConfiguration(
            direction=AxisEnum.z,
            size=vacuum,
        )

        # Get direction from kwargs or use default
        direction = kwargs.pop('direction', AxisEnum.z)
        
        return cls(
            stack_components=[atomic_layers, vacuum_config],
            xy_supercell_matrix=xy_supercell_matrix,
            direction=direction,
            **kwargs,
        )
