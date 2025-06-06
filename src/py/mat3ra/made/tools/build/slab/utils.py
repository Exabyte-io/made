from typing import Union, List
import numpy as np

from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.miller_indices import (
    MillerIndicesSchema,
)

from mat3ra.made.material import Material
from .termination import Termination
from ...convert import to_pymatgen
from ...third_party import PymatgenSlab, PymatgenSlabGenerator, label_pymatgen_slab_termination


def generate_pymatgen_slabs(
    crystal: Material,
    miller_indices: Union[MillerIndicesSchema, List[int]] = (0, 0, 1),
    min_slab_size=1,
    min_vacuum_size=0,
    in_unit_planes: bool = True,
    make_primitive: bool = False,
    symmetrize: bool = False,
) -> List[PymatgenSlab]:  # type: ignore
    # Extract actual values from MillerIndicesSchema if needed
    if isinstance(miller_indices, MillerIndicesSchema):
        miller_values = list(miller_indices.root)
    else:
        miller_values = miller_indices

    generator = PymatgenSlabGenerator(
        initial_structure=to_pymatgen(crystal),
        miller_index=miller_values,
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


def generate_miller_supercell_matrix(
    crystal: Material,
    miller_indices: Union[MillerIndicesSchema, List[int]] = (0, 0, 1),
    min_slab_size=1,
    min_vacuum_size=0,
    in_unit_planes: bool = True,
    make_primitive: bool = False,
    symmetrize: bool = False,
) -> List[List[int]]:

    if isinstance(miller_indices, MillerIndicesSchema):
        miller_values = list(miller_indices.root)
    else:
        miller_values = miller_indices

    generator = PymatgenSlabGenerator(
        initial_structure=to_pymatgen(crystal),
        miller_index=miller_values,
        min_slab_size=min_slab_size,
        min_vacuum_size=min_vacuum_size,
        in_unit_planes=in_unit_planes,
        primitive=make_primitive,
    )

    supercell_matrix = generator.slab_scale_factor.tolist()
    return supercell_matrix


def get_terminations(crystal: Material, miller_indices: Union[MillerIndicesSchema, List[int]]) -> List[Termination]:
    return [
        Termination.from_string(label_pymatgen_slab_termination(slab))
        for slab in generate_pymatgen_slabs(crystal, miller_indices)
    ]


def choose_termination(terminations: List[Termination], stoichiometry: str) -> Termination:
    # choose a termination by stochiometry or symmetry provided
    for termination in terminations:
        if str(termination.chemical_elements) == stoichiometry:
            return termination
    return terminations[0] if terminations else None
