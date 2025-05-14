from mat3ra.esse.models.material import MaterialSchema

from mat3ra.esse.models.materials_category.defects.configuration import BaseDefectConfigurationSchema


class BaseDefectConfiguration(BaseDefectConfigurationSchema):
    """
    Base schema for defect configurations.

    Args:
        crystal (MaterialSchema): The material object.
    """

    crystal: MaterialSchema
