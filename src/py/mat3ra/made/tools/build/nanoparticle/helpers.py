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

    repetitions_z = analyzer.number_of_layers_to_fit_nanoparticle
    xy_supercell = analyzer.xy_supercell_matrix_to_fit_nanoparticle

    slab = create_slab(
        material,
        miller_indices=orientation_z,
        number_of_layers=repetitions_z,
        xy_supercell_matrix=xy_supercell,
        vacuum=0,
    )

    nanoparticle = filter_by_condition_on_coordinates(slab, condition.condition, use_cartesian_coordinates=True)

    return nanoparticle
