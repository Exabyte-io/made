from mat3ra.made.material import Material
from .builders import PointDefectBuilder, PointDefectBuilderParameters
from .configuration import PointDefectConfiguration


def create_defect(material: Material, configuration: PointDefectConfiguration) -> Material:
    """
    Add a defect to a material.

    Args:
        material: The material to which the defect will be added.
        configuration: The configuration of the defect to be added.

    Returns:
        The material with the defect added.
    """
    builder = PointDefectBuilder(build_parameters=PointDefectBuilderParameters())
    configuration.material = material
    material = builder.get_material(configuration)
    return material
