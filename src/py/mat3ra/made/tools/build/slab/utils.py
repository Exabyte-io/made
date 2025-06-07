from typing import Union, List

from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.miller_indices import (
    MillerIndicesSchema,
)
from mat3ra.made.material import Material
from .termination import Termination
from ...convert import to_pymatgen
from ...third_party import PymatgenSlab, PymatgenSlabGenerator, label_pymatgen_slab_termination


def create_pymatgen_slab_generator(
    crystal: Material,
    miller_indices: Union[MillerIndicesSchema, List[int]],
    min_slab_size: float,
    min_vacuum_size: float,
    in_unit_planes: bool,
    make_primitive: bool,
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
    miller_indices: Union[MillerIndicesSchema, List[int]] = (0, 0, 1),
    min_slab_size: float = 1.0,
    min_vacuum_size: float = 0.0,
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


def get_terminations(crystal: Material, miller_indices: Union[MillerIndicesSchema, List[int]]) -> List[Termination]:
    DEFAULT_THICKNESS = 3
    DEFAULT_VACUUM_SIZE = 1
    slabs = generate_pymatgen_slabs(
        crystal, miller_indices, min_slab_size=DEFAULT_THICKNESS, min_vacuum_size=DEFAULT_VACUUM_SIZE, symmetrize=True
    )
    return [Termination.from_string(label_pymatgen_slab_termination(slab)) for slab in slabs]


def select_termination(terminations: List[Termination], stoichiometry: str) -> Termination:
    for termination in terminations:
        if str(termination.chemical_elements) == stoichiometry:
            return termination
    print("No stoichiometry found. Selecting the first termination.")
    return terminations[0] if terminations else None
