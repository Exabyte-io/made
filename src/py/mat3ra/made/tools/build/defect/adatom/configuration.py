from mat3ra.esse.models.core.abstract.coordinate_3d import Coordinate3dSchema
from mat3ra.esse.models.materials_category.defective_structures.two_dimensional.adatom.configuration import (
    AdatomPointDefectSchema,
    AxisEnum,
)
from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import ChemicalElement

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.slab.configuration import SlabStackConfiguration


class AdatomDefectConfiguration(SlabStackConfiguration):
    """
    Configuration for creating an adatom defect on a slab surface.

    Args:
        stack_components: List containing [slab, isolated_defect, vacuum].
    """

    type: str = "AdatomDefectConfiguration"
    direction: AxisEnum = AxisEnum.z
