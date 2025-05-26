from typing import List, Union, Optional
from pydantic import BaseModel, confloat, Field

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


# termination = CrystalLatticePlanesSchema.get_terminations()


class AtomicLayers(AtomicLayersSchema, CrystalLatticePlanes):
    # terminations: List[Termination]
    pass


class AtomicLayersUnique(AtomicLayersUniqueSchema, CrystalLatticePlanes):
    # terminations: List[Termination]
    pass


class AtomicLayersUniqueRepeated(AtomicLayersUniqueRepeatedSchema, CrystalLatticePlanes):
    crystal: Material
    termination_top: Termination

    def _generate_pymatgen_slabs(self, **kwargs) -> List[PymatgenSlab]:
        return super()._generate_pymatgen_slabs(
            min_slab_size=self.number_of_repetitions,
            min_vacuum_size=kwargs.get("min_vacuum_size", 1),
            in_unit_planes=kwargs.get("in_unit_planes", True),
            make_primitive=kwargs.get("make_primitive", False),
            symmetrize=kwargs.get("symmetrize", False),
        )


# TODO: fix pydantic generation in ESSE and use
class VacuumConfiguration(BaseModel):
    """
    Configuration for adding vacuum to a material.

    Args:
        direction (AxisEnum): The direction to add vacuum along (x, y, or z).
        size (float): Size of the vacuum gap in angstroms. Defaults to 10.0.
        is_orthogonal (bool): Whether the vacuum slab is orthogonal to the specified direction. Defaults to True.
    """

    direction: AxisEnum
    size: float = 10.0
    is_orthogonal: bool = True


class SlabConfiguration(SlabSchema, BaseConfigurationPydantic):
    type: str = "SlabConfiguration"
    stack_components: List[Union[AtomicLayersUniqueRepeated, VacuumConfiguration]]
    supercell_xy: List[List[int]]
    direction: AxisEnum = AxisEnum.z
