from typing import List, Union

from mat3ra.esse.models.core.abstract.vector_3d import Vector3dSchema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.miller_indices import (
    MillerIndicesSchema,
)
from mat3ra.made.material import Material
from pydantic import BaseModel

from ..convert import from_pymatgen, to_pymatgen
from ..third_party import PymatgenSlab, PymatgenSlabGenerator, label_pymatgen_slab_termination
from .lattice import LatticeMaterialAnalyzer
from .termination import Termination


def select_slab_with_termination_by_formula(slabs: List[PymatgenSlab], termination: Termination) -> PymatgenSlab:
    for slab in slabs:
        if Termination.from_string(label_pymatgen_slab_termination(slab)).formula == termination.formula:
            return slab
    raise ValueError(f"No slab found with termination: {termination}")


class TerminationHolder(BaseModel):
    termination_with_vacuum: Termination
    termination_without_vacuum: Termination
    shift: float


class CrystalLatticePlanesMaterialAnalyzer(LatticeMaterialAnalyzer):
    DEFAULT_THICKNESS = 3
    DEFAULT_VACUUM_SIZE = 1
    DEFAULT_SYMMETRIZE = False

    def __init__(self, material: Material, miller_indices: Union[MillerIndicesSchema, List[int]]):
        super().__init__(material)
        self.miller_indices = miller_indices
        self.termination_holder = self.get_termination_holders()

    def get_pymatgen_slab_generator(
        self,
        min_slab_size: float = DEFAULT_VACUUM_SIZE,
        min_vacuum_size: float = DEFAULT_VACUUM_SIZE,
        in_unit_planes: bool = True,
        make_primitive: bool = False,
    ):
        return PymatgenSlabGenerator(
            initial_structure=to_pymatgen(self.material_with_conventional_lattice),
            miller_index=self.miller_indices,
            min_slab_size=min_slab_size,
            min_vacuum_size=min_vacuum_size,
            in_unit_planes=in_unit_planes,
            primitive=make_primitive,
        )

    @property
    def pymatgen_slab_generator_with_vacuum(self) -> PymatgenSlabGenerator:
        return self.get_pymatgen_slab_generator()

    @property
    def pymatgen_slab_generator_without_vacuum(self) -> PymatgenSlabGenerator:
        return self.get_pymatgen_slab_generator(min_vacuum_size=0)

    @property
    def all_planes_as_pymatgen_slabs_with_vacuum(self) -> List[PymatgenSlab]:
        return self.pymatgen_slab_generator_with_vacuum.get_slabs(symmetrize=self.DEFAULT_SYMMETRIZE)

    @property
    def all_planes_as_pymatgen_slabs_without_vacuum(self) -> List[PymatgenSlab]:
        return self.pymatgen_slab_generator_without_vacuum.get_slabs(symmetrize=self.DEFAULT_SYMMETRIZE)

    def get_termination_holders(self):
        termination_holders = []
        slabs_with_vacuum = self.all_planes_as_pymatgen_slabs_with_vacuum
        slabs_without_vacuum = self.all_planes_as_pymatgen_slabs_without_vacuum

        for slab_with_vacuum in slabs_with_vacuum:
            termination_with_vacuum_string = label_pymatgen_slab_termination(slab_with_vacuum)
            termination_with_vacuum = Termination.from_string(termination_with_vacuum_string)

            matching_slab_without_vacuum = select_slab_with_termination_by_formula(
                slabs_without_vacuum, termination_with_vacuum
            )

            termination_without_vacuum_string = label_pymatgen_slab_termination(matching_slab_without_vacuum)
            termination_without_vacuum = Termination.from_string(termination_without_vacuum_string)

            termination_holders.append(
                TerminationHolder(
                    termination_with_vacuum=termination_with_vacuum,
                    termination_without_vacuum=termination_without_vacuum,
                    shift=slab_with_vacuum.shift,
                )
            )

        return termination_holders

    @property
    def terminations_with_vacuum(self) -> List[Termination]:
        return [holder.termination_with_vacuum for holder in self.termination_holders]

    @property
    def terminations_without_vacuum(self) -> List[Termination]:
        return [holder.termination_without_vacuum for holder in self.termination_holders]

    @property
    def terminations(self):
        return self.terminations_with_vacuum

    @property
    def default_termination(self) -> Termination:
        return self.terminations[0]

    @property
    def miller_supercell_matrix(self) -> Matrix3D:  # TODO: import from mat3ra.esse.models.core
        return self.pymatgen_slab_generator_without_vacuum.slab_scale_factor.tolist()

    def get_material_with_termination_without_vacuum(self, termination: Termination) -> Material:
        # Note: termination is passed with vacuum, but we return the material without vacuum
        holder = next((h for h in self.termination_holders if h.termination_with_vacuum == termination), None)
        if holder is None:
            raise ValueError(f"Termination {termination} not found.")
        slab = self.all_planes_as_pymatgen_slabs_without_vacuum[holder.index]
        return Material.create(from_pymatgen(slab))

    def get_translation_vector_for_termination_without_vacuum(self, termination: Termination) -> Vector3dSchema:
        holder = next((h for h in self.termination_holders if h.termination_with_vacuum == termination), None)
        if holder is None:
            raise ValueError(f"Termination {termination} not found.")
        return holder.shift
