from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import AtomSchema

from mat3ra.made.material import Material
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.point.builders import AtomAtCoordinateConfiguration
from mat3ra.made.tools.build.passivation.analyzer import SurfacePassivationMaterialAnalyzer
from mat3ra.made.tools.build.passivation.builders import PassivationBuilder
from mat3ra.made.tools.build.passivation.configuration import PassivationConfiguration


def create_passivated_surface(
    material: MaterialWithBuildMetadata,
    passivant: str = "H",
    bond_length: float = 1.0,
    shadowing_radius: float = 2.5,
    depth: float = 5.0,
) -> Material:
    """
    Create a passivated material by adding passivant atoms to the surface of the provided material.

    Args:
        material (Material): The material to be passivated.
        passivant (str): The chemical symbol of the passivating atom (default is 'H').
        bond_length (float): The bond length for the passivation (default is 1.0).

    Returns:
        Material: The passivated material.
    """
    analyzer = SurfacePassivationMaterialAnalyzer(
        material=material, passivant=passivant, bond_length=bond_length, shadowing_radius=shadowing_radius, depth=depth
    )

    passivant_coordinates = analyzer.get_passivant_coordinates()
    passivant_configs = [
        AtomAtCoordinateConfiguration(
            crystal=material, element=AtomSchema(chemical_element=passivant), coordinate=coord
        )
        for coord in passivant_coordinates
    ]
    config = PassivationConfiguration(
        merge_components=[material, *passivant_configs], passivant=passivant, bond_length=bond_length
    )

    builder = PassivationBuilder()
    return builder.get_material(config)
