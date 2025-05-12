from pydantic import BaseModel

from mat3ra.esse.models.materials_category.single_material.material import MaterialSchema


class BaseDefectConfigurationSchema(BaseModel):
    """
    Base schema for defect configurations.

    Args:
        crystal (MaterialSchema): The material object.
    """
    crystal: MaterialSchema 