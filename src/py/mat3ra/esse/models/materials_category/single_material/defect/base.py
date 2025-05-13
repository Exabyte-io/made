from mat3ra.esse.models.materials_category.single_material.material import MaterialSchema
from pydantic import BaseModel


class BaseDefectConfigurationSchema(BaseModel):
    """
    Base schema for defect configurations.

    Args:
        crystal (MaterialSchema): The material object.
    """

    crystal: MaterialSchema
