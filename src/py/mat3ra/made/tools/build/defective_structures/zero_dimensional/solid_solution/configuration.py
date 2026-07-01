from typing import List, Optional, Union

from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.solid_solution.configuration import (
    SolidSolutionConfigurationSchema,
)
from mat3ra.made.material import Material

from .....build_components import BaseConfigurationPydantic, MaterialWithBuildMetadata


class SolidSolutionConfiguration(BaseConfigurationPydantic, SolidSolutionConfigurationSchema):
    type: str = "SolidSolutionConfiguration"
    crystal: Union[Material, MaterialWithBuildMetadata]
    source_element: str
    target_element: str
    concentration: float
    supercell_dimensions: List[int] = [1, 1, 1]
    seed: Optional[int] = None
    site_selection_method: str = "random"

