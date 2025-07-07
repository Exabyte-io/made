from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer

from mat3ra.made.tools.build.defect.configuration import SlabDefectConfigurationLegacy
from mat3ra.made.tools.build.defect.slab.builders import SlabDefectBuilder
from mat3ra.made.tools.build.defect.slab.configuration import SlabDefectConfiguration
from mat3ra.made.tools.build.slab.builders import SlabWithAdditionalLayersBuilder


def create_slab_defect(slab: Material, isolated_defect: Material) -> Material:
    """

    Create a slab defect by merging an isolated defect with a slab material.
    Args:
        slab: The original slab material.
    Returns:
        Material: The new slab material with additional layers and vacuum if needed.
    """
    analyzer = SlabMaterialAnalyzer(material=slab, auto_add_vacuum=True)
    slab_config = analyzer.slab_with_additional_layers_config()
    adjusted_slab = analyzer.slab_with_added_vacuum()

    builder = SlabDefectBuilder()
    config = SlabDefectConfiguration.from_paramteres(
        slab=adjusted_slab,
        isolated_defect=isolated_defect,
    )

    return builder.get_material(config)
