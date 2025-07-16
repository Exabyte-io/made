from typing import Union, List

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.other import get_closest_site_id_from_coordinate, get_chemical_formula
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.island.helpers import CoordinateConditionType
from mat3ra.made.tools.build.nanoparticle.analyzer import NanoparticleMaterialAnalyzer
from mat3ra.made.tools.build.nanoparticle.builders import ASEBasedNanoparticleBuilder
from mat3ra.made.tools.build.nanoparticle.configuration import ASEBasedNanoparticleConfiguration
from mat3ra.made.tools.build.nanoparticle.enums import NanoparticleShapesEnum
from mat3ra.made.tools.build.slab.helpers import create_slab
from mat3ra.made.tools.modify import filter_by_condition_on_coordinates
from mat3ra.made.tools.utils.coordinate import SphereCoordinateCondition


def create_nanoparticle_from_material(
    material: Material,
    condition: Union[
        SphereCoordinateCondition,
        CoordinateConditionType,
    ] = SphereCoordinateCondition(),
    orientation_z: List[int] = [0, 0, 1],
    center_around_atom: bool = True,
    use_cartesian_coordinates: bool = True,
    vacuum_padding: float = 10.0,
) -> MaterialWithBuildMetadata:
    """
    Create a nanoparticle by applying a coordinate condition to a crystal.

    Args:
        material: The crystal material from which the nanoparticle is created.
        condition: The coordinate condition that defines the nanoparticle shape.
        orientation_z: The orientation of the crystallographic axis in the z-direction.
        center_around_atom: Whether to center the condition around an atom close to the center in the material.
        use_cartesian_coordinates: Whether the condition supplied is in Cartesian coordinates.
        vacuum_padding: Padding for vacuum space around the nanoparticle.

    Returns:
        MaterialWithBuildMetadata: The nanoparticle material with build metadata.
    """

    analyzer = NanoparticleMaterialAnalyzer(
        material=material,
        orientation_z=orientation_z,
        vacuum_padding=vacuum_padding,
        coordinate_condition=condition,
    )

    repetitions_z = analyzer.number_of_layers_to_fit_nanoparticle
    xy_supercell = analyzer.xy_supercell_matrix_to_fit_nanoparticle

    slab = create_slab(
        material,
        miller_indices=orientation_z,
        number_of_layers=repetitions_z,
        xy_supercell_matrix=xy_supercell,
        vacuum=0,
    )

    # Find an atom close to center, to apply the condition around it
    if center_around_atom:
        slab.to_cartesian()
        center_coordinate = slab.basis.cell.convert_point_to_cartesian([0.5, 0.5, 0.5])
        center_id_at_site = get_closest_site_id_from_coordinate(slab, center_coordinate, use_cartesian_coordinates=True)
        center_coordinate = slab.basis.coordinates.get_element_value_by_index(center_id_at_site)
        condition.center_coordinate = center_coordinate

    nanoparticle = filter_by_condition_on_coordinates(
        slab, condition.condition, use_cartesian_coordinates=use_cartesian_coordinates
    )

    formula = get_chemical_formula(nanoparticle)
    coordinate_condition_name = condition.__class__.__name__.replace("CoordinateCondition", "")

    nanoparticle.name = f"{formula} {coordinate_condition_name} Nanoparticle"

    return nanoparticle


def create_nanoparticle_by_shape(
    crystal: Material,
    shape: NanoparticleShapesEnum,
    parameters: dict,
):
    """
    Create a nanoparticle from a crystal material by specifying its shape and parameters.

    Provide the crystal which first element nad lattice constant will be used to create the nanoparticle.

    """

    # TODO: should done in analyzer
    element = crystal.basis.elements.values[0]
    lattice_constant = crystal.lattice.a

    return create_nanoparticle_by_shape_from_element(
        element=element, shape=shape, lattice_constant=lattice_constant, parameters=parameters
    )


def create_nanoparticle_by_shape_from_element(
    element: str,
    lattice_constant: float,
    shape: NanoparticleShapesEnum,
    parameters: dict,
):
    """
    Create a nanoparticle from a specified element, lattice constant, shape, and parameters.
    The shape is defined by the NanoparticleShapesEnum, and parameters are passed directly to the ASE constructor.
    """
    if "latticeconstant" not in parameters:
        parameters["latticeconstant"] = lattice_constant

    config = ASEBasedNanoparticleConfiguration(shape=shape, parameters=parameters, element=element)
    builder = ASEBasedNanoparticleBuilder()
    nanoparticles = builder.get_materials(config)
    return nanoparticles[0]
