from typing import List, Union

from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.miller_indices import (
    MillerIndicesSchema,
)
from pydantic import BaseModel

from mat3ra.made.material import Material
from src.py.mat3ra.esse.models.core.abstract.vector_3d import Vector3dSchema
from .lattice import LatticeMaterialAnalyzer
from .termination import Termination
from ..convert import from_pymatgen, to_pymatgen
from ..third_party import PymatgenSlab, label_pymatgen_slab_termination, PymatgenSlabGenerator


def generate_miller_supercell_matrix(
    crystal: Material,
    miller_indices: Union[MillerIndicesSchema, List[int]] = (0, 0, 1),
    min_slab_size: float = 1.0,
    min_vacuum_size: float = 0.0,
    in_unit_planes: bool = True,
    make_primitive: bool = False,
) -> List[List[int]]:
    generator = create_pymatgen_slab_generator(
        crystal, miller_indices, min_slab_size, min_vacuum_size, in_unit_planes, make_primitive
    )
    return generator.slab_scale_factor.tolist()


def create_pymatgen_slab_generator(
    crystal: Material,
    miller_indices: Union[MillerIndicesSchema, List[int]],
    min_slab_size: float = 1.0,
    min_vacuum_size: float = 0.0,
    in_unit_planes: bool = True,
    make_primitive: bool = False,
) -> PymatgenSlabGenerator:
    if isinstance(miller_indices, MillerIndicesSchema):
        miller_values = list(miller_indices.root)
    else:
        miller_values = miller_indices
    return PymatgenSlabGenerator(
        initial_structure=to_pymatgen(crystal),
        miller_index=miller_values,
        min_slab_size=min_slab_size,
        min_vacuum_size=min_vacuum_size,
        in_unit_planes=in_unit_planes,
        primitive=make_primitive,
    )


def generate_pymatgen_slabs(
    crystal: Material,
    miller_indices: Union[MillerIndicesSchema, List[int]],
    min_slab_size: float,
    min_vacuum_size: float,
    in_unit_planes: bool = True,
    make_primitive: bool = False,
    symmetrize: bool = False,
) -> List[PymatgenSlab]:  # type: ignore
    generator = create_pymatgen_slab_generator(
        crystal, miller_indices, min_slab_size, min_vacuum_size, in_unit_planes, make_primitive
    )
    return generator.get_slabs(
        symmetrize=symmetrize,
    )


class TerminationHolder(BaseModel):
    termination_with_vacuum: Termination
    termination_without_vacuum: Termination
    shift: float


class CrystalLatticePlanesMaterialAnalyzer(LatticeMaterialAnalyzer):
    DEFAULT_THICKNESS = 3
    DEFAULT_VACUUM_SIZE = 1

    def __init__(self, material: Material, miller_indices: Union[MillerIndicesSchema, List[int]]):
        super().__init__(material)
        self.miller_indices = miller_indices
        self.init_termination_holders()

    def pymatgen_slab_generator_with(self):
        pass

    @property
    def all_planes_as_pymatgen_slabs_with_vacuum(self) -> List[PymatgenSlab]:
        # Implement here
        return generate_pymatgen_slabs(
            crystal=self.material_with_conventional_lattice,
            miller_indices=self.miller_indices,
            min_slab_size=self.DEFAULT_THICKNESS,
            min_vacuum_size=self.DEFAULT_VACUUM_SIZE,
        )

    @property
    def all_planes_as_pymatgen_slabs_without_vacuum(self) -> List[PymatgenSlab]:
        # Implement here
        return generate_pymatgen_slabs(
            crystal=self.material_with_conventional_lattice,
            miller_indices=self.miller_indices,
            min_slab_size=1,
            min_vacuum_size=0,
        )

    def get_all_planes_as_pymatgen_slabs_symmetrized(self) -> List[PymatgenSlab]:
        pass

    def init_termination_holders(self):
        termination_holders = []
        slabs_with_vacuum = self.all_planes_as_pymatgen_slabs_with_vacuum
        slabs_without_vacuum = self.all_planes_as_pymatgen_slabs_without_vacuum

        for slab_with_vacuum in slabs_with_vacuum:
            termination_with_vacuum = Termination.from_string(label_pymatgen_slab_termination(slab_with_vacuum))

            matching_slab_without_vacuum = next(
                (
                    slab
                    for slab in slabs_without_vacuum
                    if Termination.from_string(label_pymatgen_slab_termination(slab)).formula
                    == termination_with_vacuum.formula
                ),
                None,
            )
            termination_without_vacuum = (
                Termination.from_string(label_pymatgen_slab_termination(matching_slab_without_vacuum))
                if matching_slab_without_vacuum
                else None
            )

            if matching_slab_without_vacuum is None:
                raise ValueError(f"No matching non-vacuum slab found for termination: {termination_with_vacuum}")

            termination_holders.append(
                TerminationHolder(
                    termination_with_vacuum=termination_with_vacuum,
                    termination_without_vacuum=termination_without_vacuum,
                    shift=slab_with_vacuum.shift,
                )
            )

        self.termination_holders = termination_holders  # Save for reuse
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
    def miller_supercell(self) -> List[List[int]]:
        return generate_miller_supercell_matrix(self.material, self.miller_indices)

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
