from typing import List, Union

from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.atomic_layers import AtomicLayersSchema
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.atomic_layers_unique import (
    AtomicLayersUniqueSchema,
)
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.atomic_layers_unique_repeated import (
    AtomicLayersUniqueRepeatedSchema,
)
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import (
    SlabSchema,
    VacuumSchema,
    AxisEnum,
)
from ...third_party import PymatgenSlabGenerator, label_pymatgen_slab_termination
from mat3ra.made.material import Material
from .termination import Termination
from .. import BaseConfigurationPydantic
from ...convert import to_pymatgen, from_pymatgen
from ...third_party import PymatgenSlab, PymatgenSpacegroupAnalyzer


class CrystalLatticePlanes(CrystalLatticePlanesSchema):

    crystal: Material
    _pymatgen_slabs: List[PymatgenSlab]

    def __init__(self, **kwargs):
        crystal = kwargs.pop("crystal", None) or Material.create_default()
        use_conventional = kwargs.get(
            "use_conventional_cell", CrystalLatticePlanesSchema.model_fields["use_conventional_cell"].default
        )

        # TODO: use LatticeAnalyzer to get the conventional cell
        pymatgen_structure = to_pymatgen(crystal)
        if use_conventional:
            pymatgen_structure = PymatgenSpacegroupAnalyzer(pymatgen_structure).get_conventional_standard_structure()

        crystal_config = from_pymatgen(pymatgen_structure)
        kwargs["crystal"] = Material.create(crystal_config)
        super().__init__(**kwargs)
        self._pymatgen_slabs = self._get_pymatgen_slabs(crystal)

    # TODO: add specific parameters
    def _generate_pymatgen_slabs(self, symmetrize) -> List[PymatgenSlab]:  # type: ignore
        generator = PymatgenSlabGenerator(
            initial_structure=to_pymatgen(self.crystal),
            miller_index=self.miller_indices,
            min_slab_size=1,  # TODO: change
            min_vacuum_size=1,  # TODO: change
            in_unit_planes=True,
            primitive=False,
        )
        raw_slabs = generator.get_slabs(
            # We need to preserve symmetric slabs for different terminations at the surface
            symmetrize=symmetrize
        )

        return raw_slabs

    def get_terminations(self):
        """
        Get the terminations for the slab.
        """
        return [Termination.from_string(label_pymatgen_slab_termination(slab)) for slab in self._pymatgen_slabs]


# termination = CrystalLatticePlanesSchema.get_terminations()


class AtomicLayers(AtomicLayersSchema, CrystalLatticePlanes):
    # terminations: List[Termination]
    pass


class AtomicLayersUnique(AtomicLayersUniqueSchema, CrystalLatticePlanes):
    # terminations: List[Termination]
    pass


class AtomicLayersUniqueRepeated(AtomicLayersUniqueRepeatedSchema, CrystalLatticePlanes):
    termination_top: Termination


class VacuumConfiguration(VacuumSchema):
    pass


class SlabConfiguration(SlabSchema, BaseConfigurationPydantic):
    type: str = "SlabConfiguration"
    stack_components: List[Union[AtomicLayersUniqueRepeated, VacuumConfiguration]]
    supercell_matrix: List[List[int]]
    direction: AxisEnum = AxisEnum.z
