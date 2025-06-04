from mat3ra.esse.models.material.primitive.two_dimensional.vacuum import VacuumConfigurationSchema

from mat3ra.made.material import Material


class VacuumConfiguration(VacuumConfigurationSchema):
    size: float
    crystal: Material
