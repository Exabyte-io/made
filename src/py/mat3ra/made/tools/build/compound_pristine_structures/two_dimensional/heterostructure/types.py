from typing import List, Optional, Tuple, Union

from mat3ra.code.entity import InMemoryEntityPydantic
from pydantic import Field

from mat3ra.made.material import Material
from .....build_components import MaterialWithBuildMetadata


class StackComponentDict(InMemoryEntityPydantic):
    """
    Pydantic model for heterostructure stack component configurations.

    Required fields:
        crystal: The crystal material to create a slab from
        miller_indices: Miller indices for the slab surface as (h, k, l)
        thickness: Number of layers in the slab

    Optional fields:
        xy_supercell_matrix: Optional supercell matrix for the xy plane
    """

    crystal: Union[Material, MaterialWithBuildMetadata] = Field(..., description="Crystal material for the slab")
    miller_indices: Tuple[int, int, int] = Field(..., description="Miller indices for the slab surface")
    thickness: int = Field(..., gt=0, description="Number of layers in the slab")
    xy_supercell_matrix: Optional[List[List[int]]] = Field(None, description="Optional xy supercell matrix")
