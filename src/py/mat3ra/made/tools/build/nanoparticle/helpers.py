from typing import Union, List

from mat3ra.made.material import Material
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.island.helpers import CoordinateConditionType
from mat3ra.made.tools.build.nanoparticle.analyzer import NanoparticleMaterialAnalyzer
from mat3ra.made.tools.build.nanoparticle.configuration import NanoparticleConfiguration
from mat3ra.made.tools.build.slab.helpers import create_slab
from mat3ra.made.tools.build.supercell import create_supercell
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
    )

    # Create a supercell for the slab to cut nanoparticle
    supercell_to_fit_nanoparticle = analyzer.supercell_to_fit_nanpoarticle
    repetitions_z = supercell_to_fit_nanoparticle[2][2]
    xy_supercell = supercell_to_fit_nanoparticle[:2][0:2]

    # create supercell that includes that + vacuum padding on all sides
    supercell_for_vacuum = analyzer.supercell_for_vacuum

    slab = create_slab(
        material,
        miller_indices=orientation_z,
        number_of_layers=repetitions_z,
        xy_supercell_matrix=xy_supercell,
    )

    config = NanoparticleConfiguration(merge_components=slab)
