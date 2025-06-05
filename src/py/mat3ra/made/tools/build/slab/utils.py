from typing import Union, List

from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.miller_indices import MillerIndicesSchema

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
    # Extract actual values from MillerIndicesSchema if needed
    print(f"DEBUG: miller_indices type: {type(miller_indices)}")
    print(f"DEBUG: isinstance MillerIndicesSchema check: {isinstance(miller_indices, MillerIndicesSchema)}")
    print(f"DEBUG: miller_indices value: {miller_indices}")
    
    if isinstance(miller_indices, MillerIndicesSchema):
        miller_values = list(miller_indices.root)
        print(f"DEBUG: Extracted miller_values from MillerIndicesSchema: {miller_values}, type: {type(miller_values)}")
    else:
        miller_values = miller_indices
        print(f"DEBUG: Using miller_indices directly: {miller_values}, type: {type(miller_values)}")

    print(f"DEBUG: About to pass miller_index={miller_values} to PymatgenSlabGenerator")

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


def calculate_rotation_matrix(crystal, miller_supercell_material):
    # Implement logic to calculate the rotation matrix based on crystal and miller_supercell_material
    return [[1, 0, 0], [0, 1, 0], [0, 0, 1]]  # Identity matrix as a placeholder


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
