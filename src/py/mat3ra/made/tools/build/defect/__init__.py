from enum import Enum
from typing import Optional, List, Literal, Union

from pydantic import BaseModel

from src.py.mat3ra.made.material import Material
from ...build import BaseBuilder


class PointDefectTypeEnum(str, Enum):
    VACANCY = "vacancy"
    SUBSTITUTION = "substitution"
    INTERSTITIAL = "interstitial"


class BaseDefectConfiguration(BaseModel):
    material: Material
    # defect_type type can be an Enum for a specific defect class (for point defect, 2d defect, etc.)
    defect_type: Union[PointDefectTypeEnum, None] = None


class PointDefectConfiguration(BaseDefectConfiguration):
    defect_type: PointDefectTypeEnum = PointDefectTypeEnum.VACANCY
    # TODO: should come from Enum of elements
    specie: Optional[str] = None
    min_distance: Optional[float] = 0.0


class PointDefectBuilderParameters(BaseModel):
    target_site: int = 0
    # TODO: import coordinate type from ESSE
    position: Optional[List[int]] = None
    center_defect: bool = False
    resize_cell_matrix: Union[List[int], Literal["auto"]] = "auto"


class PointDefectBuilder(BaseBuilder):
    """
    Builder class for generating point defects.
    """

    _BuildParametersType = PointDefectBuilderParameters
    _DefaultBuildParameters = PointDefectBuilderParameters()
    _GeneratedItemType = ...
    _ConfigurationType = PointDefectConfiguration
    _SelectorParametersType = None
    _PostProcessParametersType = None