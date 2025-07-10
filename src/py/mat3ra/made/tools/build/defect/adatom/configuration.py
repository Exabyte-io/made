from typing import List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.defective_structures.two_dimensional.adatom.configuration import (
    AdatomPointDefectSchema,
)
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.slab.configuration import SlabStackConfiguration
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration


class AdatomDefectConfiguration(MergeConfiguration, AdatomPointDefectSchema):
    """
    Configuration for creating an adatom defect on a slab surface.S

    Args:
        merge_components: List containing [slab, isolated_defect].
    """

    type: str = "AdatomDefectConfiguration"


class CrystalSiteAdatomConfiguration(SlabStackConfiguration):
    """
    Configuration for creating an adatom defect using crystal site placement with SlabStackBuilder.

    This configuration uses the new approach where:
    1. Original slab (without vacuum)
    2. Added component (single atom in lattice created using recreate_slab_with_fractional_layers)
    3. Vacuum layer

    Args:
        stack_components: List containing [slab_without_vacuum, added_component, vacuum].
        coordinate: Crystal site coordinate for the adatom.
        element: Chemical element of the adatom.
    """

    type: str = "CrystalSiteAdatomConfiguration"
    coordinate: List[float]
    element: str

    @classmethod
    def from_parameters(cls, slab: Material, coordinate: List[float], element: str) -> "CrystalSiteAdatomConfiguration":
        """
        Create a CrystalSiteAdatomConfiguration from slab, coordinate and element.

        Args:
            slab: The original slab material
            coordinate: Crystal coordinate for the adatom
            element: Chemical element of the adatom

        Returns:
            CrystalSiteAdatomConfiguration: The created configuration
        """

        # Get slab without vacuum
        analyzer = SlabMaterialAnalyzer(material=slab)
        slab_without_vacuum = analyzer.get_slab_configuration_with_no_vacuum()

        # Get vacuum configuration
        vacuum_config = analyzer.get_slab_vacuum_configuration()

        # Initialize with placeholder for added_component (will be filled by builder)
        stack_components = [slab_without_vacuum, None, vacuum_config]

        return cls(stack_components=stack_components, direction=AxisEnum.z, coordinate=coordinate, element=element)
