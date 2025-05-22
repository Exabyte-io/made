from typing import List, Union

from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.atomic_layers import AtomicLayersSchema
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import (
    SlabSchema,
    VacuumSchema,
    AxisEnum,
)


from mat3ra.made.material import Material
from .termination import Termination
from .. import BaseConfigurationPydantic
from ...convert import to_pymatgen, from_pymatgen
from ...third_party import PymatgenSpacegroupAnalyzer


class CrystalLatticePlanes(CrystalLatticePlanesSchema):
    crystal: Material
    miller_indices: List[int] = [1, 1, 1]


# termination = CrystalLatticePlanesSchema.termination


class AtomicLayers(AtomicLayersSchema, CrystalLatticePlanes):
    # terminations: List[Termination]
    pass


class VacuumConfiguration(VacuumSchema):
    direction: float
    size: float


class SlabConfiguration(SlabSchema, BaseConfigurationPydantic):
    """
    Configuration for building a slab.

    Args:
        bulk (Material): The bulk material.
        miller_indices (Tuple[int, int, int]): The Miller indices of the slab.
        thickness (int): The thickness of the slab.
        vacuum (float): The vacuum thickness, in Angstroms.
        xy_supercell_matrix (List[List[int]]): The supercell matrix for the xy plane.
        use_conventional_cell (bool): Whether to use the conventional cell.
        use_orthogonal_z (bool): Whether to use orthogonal z.
        make_primitive (bool): Whether to try to find primitive cell for the created slab.
    """

    type: str = "SlabConfiguration"
    stack_components: List[Union[AtomicLayers, VacuumConfiguration]]
    supercell_matrix: List[List[int]]
    direction: AxisEnum = "z"

    def __init__(self, **kwargs):
        bulk = kwargs.pop("bulk", None) or Material.create_default()
        use_conventional = kwargs.get("use_conventional_cell", True)

        bulk_pmg = to_pymatgen(bulk)
        if use_conventional:
            bulk_pmg = PymatgenSpacegroupAnalyzer(bulk_pmg).get_conventional_standard_structure()

        bulk_config = from_pymatgen(bulk_pmg)
        kwargs["bulk"] = Material.create(bulk_config)

        super().__init__(**kwargs)
