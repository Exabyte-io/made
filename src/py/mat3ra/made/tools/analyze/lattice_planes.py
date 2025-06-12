from typing import List, Tuple, Union

from mat3ra.code.vector import Vector3D
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.entities import Termination, TerminationHolder

from ..build.slab.entities import MillerIndices
from ..convert import from_pymatgen, to_pymatgen
from ..third_party import PymatgenSlab, PymatgenSlabGenerator, label_pymatgen_slab_termination
from .lattice import LatticeMaterialAnalyzer


def select_slab_with_termination_by_formula(slabs: List[PymatgenSlab], termination: Termination) -> PymatgenSlab:
    for slab in slabs:
        if Termination.from_string(label_pymatgen_slab_termination(slab)).formula == termination.formula:
            return slab
    raise ValueError(f"No slab found with termination: {termination}")


class CrystalLatticePlanesMaterialAnalyzer(LatticeMaterialAnalyzer):
    DEFAULT_THICKNESS: int = 3
    DEFAULT_VACUUM_SIZE: int = 1
    DEFAULT_SYMMETRIZE: bool = False
    miller_indices: Union[List[int], Tuple[int, int, int]]

    def get_pymatgen_slab_generator(
        self,
        min_slab_size: float = DEFAULT_THICKNESS,
        min_vacuum_size: float = DEFAULT_VACUUM_SIZE,
        in_unit_planes: bool = True,
        make_primitive: bool = False,
    ):
        return PymatgenSlabGenerator(
            initial_structure=to_pymatgen(self.material_with_conventional_lattice),
            miller_index=MillerIndices(root=self.miller_indices).to_tuple(),
            min_slab_size=min_slab_size,
            min_vacuum_size=min_vacuum_size,
            in_unit_planes=in_unit_planes,
            primitive=make_primitive,
        )

    @property
    def pymatgen_slab_generator_with_vacuum(self) -> PymatgenSlabGenerator:
        return self.get_pymatgen_slab_generator(min_slab_size=self.DEFAULT_THICKNESS, min_vacuum_size=1)

    @property
    def pymatgen_slab_generator_without_vacuum(self) -> PymatgenSlabGenerator:
        return self.get_pymatgen_slab_generator(min_slab_size=1, min_vacuum_size=0)

    @property
    def all_planes_as_pymatgen_slabs_with_vacuum(self) -> List[PymatgenSlab]:
        return [
            *self.pymatgen_slab_generator_with_vacuum.get_slabs(symmetrize=True),
            *self.pymatgen_slab_generator_with_vacuum.get_slabs(symmetrize=False),
        ]

    @property
    def all_planes_as_pymatgen_slabs_without_vacuum(self) -> List[PymatgenSlab]:
        return [
            *self.pymatgen_slab_generator_without_vacuum.get_slabs(symmetrize=self.DEFAULT_SYMMETRIZE),
            *self.pymatgen_slab_generator_without_vacuum.get_slabs(symmetrize=not self.DEFAULT_SYMMETRIZE),
        ]

    @property
    def termination_holders(self):
        termination_holders = []
        slabs_with_vacuum = self.all_planes_as_pymatgen_slabs_with_vacuum
        slabs_without_vacuum = self.all_planes_as_pymatgen_slabs_without_vacuum

        for slab_with_vacuum in slabs_with_vacuum:
            termination_with_vacuum_string = label_pymatgen_slab_termination(slab_with_vacuum)
            termination_with_vacuum = Termination.from_string(termination_with_vacuum_string)
            termination_without_vacuum = None
            shift_without_vacuum = None
            try:
                matching_slab_without_vacuum = select_slab_with_termination_by_formula(
                    slabs_without_vacuum, termination_with_vacuum
                )

                termination_without_vacuum_string = label_pymatgen_slab_termination(matching_slab_without_vacuum)
                termination_without_vacuum = Termination.from_string(termination_without_vacuum_string)
                shift_without_vacuum = matching_slab_without_vacuum.shift
            except:
                pass
            termination_holders.append(
                TerminationHolder(
                    termination_with_vacuum=termination_with_vacuum,
                    termination_without_vacuum=termination_without_vacuum,
                    shift_with_vacuum=slab_with_vacuum.shift,
                    shift_without_vacuum=shift_without_vacuum,
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
    def miller_supercell_matrix(self) -> Matrix3x3Schema:
        return self.pymatgen_slab_generator_without_vacuum.slab_scale_factor.tolist()

    def get_material_with_termination_without_vacuum(self, termination: Termination) -> Material:
        # Note: termination is passed with vacuum, but we return the material without vacuum
        holder = next((h for h in self.termination_holders if h.termination_with_vacuum == termination), None)
        if holder is None:
            raise ValueError(f"Termination {termination} not found.")
        slab = self.all_planes_as_pymatgen_slabs_without_vacuum[
            self.terminations_without_vacuum.index(holder.termination_without_vacuum)
        ]
        return Material.create(from_pymatgen(slab))

    def get_translation_vector_for_termination_without_vacuum(self, termination: Termination) -> Vector3D:
        holder = next((h for h in self.termination_holders if h.termination_with_vacuum == termination), None)
        if holder is None:
            raise ValueError(f"Termination {termination} not found.")
        return [0.0, 0.0, holder.shift_without_vacuum]
