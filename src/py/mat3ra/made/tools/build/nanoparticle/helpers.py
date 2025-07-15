from typing import Union, List

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.island.helpers import CoordinateConditionType
from mat3ra.made.tools.build.nanoparticle.analyzer import NanoparticleMaterialAnalyzer
from mat3ra.made.tools.build.nanoparticle.configuration import NanoparticleConfiguration
from mat3ra.made.tools.build.slab.helpers import create_slab
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.made.tools.build.vacuum.builders import VacuumBuilder
from mat3ra.made.tools.modify import filter_by_condition_on_coordinates
from mat3ra.made.tools.utils.coordinate import SphereCoordinateCondition


def create_nanoparticle_by_condition(
    material: Material,
    condition: Union[
        SphereCoordinateCondition,
        CoordinateConditionType,
    ] = SphereCoordinateCondition(),
    orientation_z: List[int] = [0, 0, 1],  # direction, perpendicular to Miller indices plane
    vacuum_padding: float = 10.0,
) -> MaterialWithBuildMetadata:
    """
    Create a nanoparticle by applying a coordinate condition to a crystal.

    Args:
        material: The crystal material from which the nanoparticle is created.
        condition: The coordinate condition that defines the nanoparticle shape.
        orientation_z: The orientation of the crystallographic axis in the z-direction.
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

    # Create a supercell for the slab to cut nanoparticle+ vacuum padding on all sides
    repetitions_z = analyzer.number_of_layers_to_fit_nanoparticle
    xy_supercell = analyzer.xy_supercell_matrix_to_fit_nanoparticle

    slab = create_slab(
        material,
        miller_indices=orientation_z,
        number_of_layers=repetitions_z,
        xy_supercell_matrix=xy_supercell,
        vacuum=0,  # No vacuum in slab, we'll add it later
    )

    # # Find the center coordinate of the slab for the condition
    slab.to_cartesian()
    center_coordinate = slab.basis.cell.convert_point_to_cartesian([0.5, 0.5, 0.5])
    #
    # # Apply the coordinate condition to cut the nanoparticle
    # # Set the center of the condition to the slab center
    condition.center_position = center_coordinate
    nanoparticle = filter_by_condition_on_coordinates(slab, condition.condition, use_cartesian_coordinates=True)
    #
    # # Create a vacuum cell based on the supercell for vacuum
    # vacuum_supercell_size = supercell_for_vacuum
    # vacuum_size_x = vacuum_supercell_size[0][0] * material.lattice.a
    # vacuum_size_y = vacuum_supercell_size[1][1] * material.lattice.b
    # vacuum_size_z = vacuum_supercell_size[2][2] * material.lattice.c
    #
    # # Create a cubic vacuum cell
    # vacuum_size = max(vacuum_size_x, vacuum_size_y, vacuum_size_z)
    #
    # # Create a simple cubic lattice for vacuum
    # vacuum_material = Material.create(
    #     {
    #         "name": "Vacuum",
    #         "lattice": {
    #             "a": vacuum_size,
    #             "b": vacuum_size,
    #             "c": vacuum_size,
    #             "alpha": 90.0,
    #             "beta": 90.0,
    #             "gamma": 90.0,
    #             "units": {"length": "angstrom", "angle": "degree"},
    #             "type": "CUB",
    #         },
    #         "basis": {
    #             "elements": [],
    #             "coordinates": [],
    #             "units": "crystal",
    #         },
    #     }
    # )
    #
    # # Set the nanoparticle to have the same vacuum lattice
    # nanoparticle.set_lattice(vacuum_material.lattice)
    #
    # # Center the nanoparticle in the vacuum cell
    # nanoparticle.to_cartesian()
    # nanoparticle_coords = nanoparticle.basis.coordinates.values
    # if nanoparticle_coords:
    #     # Calculate center of mass of nanoparticle
    #     coords_array = np.array(nanoparticle_coords)
    #     center_of_mass = np.mean(coords_array, axis=0)
    #     vacuum_center = [vacuum_size / 2, vacuum_size / 2, vacuum_size / 2]
    #     translation_vector = np.array(vacuum_center) - center_of_mass
    #
    #     # Translate all coordinates
    #     translated_coords = coords_array + translation_vector
    #     nanoparticle.set_coordinates(translated_coords.tolist())
    #
    # nanoparticle.to_crystal()
    # nanoparticle.name = f"{material.name} Nanoparticle"

    return nanoparticle
