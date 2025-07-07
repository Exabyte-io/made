from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.build.slab.builders import SlabBuilder

from mat3ra.made.tools.build.defect.configuration import SlabDefectConfigurationLegacy
from mat3ra.made.tools.build.defect.slab.builders import SlabDefectBuilder
from mat3ra.made.tools.build.defect.slab.configuration import SlabDefectConfiguration


def create_slab_defect(slab: Material, isolated_defect: Material, additional_layers: int = 1) -> Material:
    """
    Create a slab defect by merging an isolated defect with a slab material.
    Args:
        slab: The original slab material.
        isolated_defect: The isolated defect material.
        additional_layers: Number of additional layers to add to the slab for defect creation.
    Returns:
        Material: The new slab material with additional layers and vacuum if needed.
    """
    analyzer = SlabMaterialAnalyzer(material=slab)
    slab_with_additional_layers_config, slab_with_original_layers_config = (
        analyzer.get_slab_with_additional_layers_configurations(
            additional_layers=additional_layers, vacuum_thickness=5.0
        )
    )

    builder = SlabBuilder()
    slab_with_additional_layers = builder.get_material(slab_with_additional_layers_config)
    slab_with_original_layers = builder.get_material(slab_with_original_layers_config)

    # defect is build with slab_with_additional_layers
    builder = SlabDefectBuilder()
    configuration = SlabDefectConfiguration(
        merge_components=[slab_with_original_layers, isolated_defect],
        merge_method=SlabDefectConfigurationLegacy.MergeMethodsEnum.ADD,
    )
    slab_with_defect = builder.get_material(configuration)
    return slab_with_defect
