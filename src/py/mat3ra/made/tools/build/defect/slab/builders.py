from pydantic import BaseModel

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.slab.configuration import SlabDefectConfiguration
from mat3ra.made.tools.build.merge import MergeBuilder

from src.py.mat3ra.made.tools.build.defect.slab.helpers import create_slab_with_additional_layers


class SlabDefectBuilderParameters(BaseModel):
    auto_add_vacuum: bool = True
    vacuum: float = 5.0


class SlabDefectBuilder(MergeBuilder):
    _ConfigurationType = SlabDefectConfiguration

    def get_slab_with_additional_layers(self, configuration: SlabDefectConfiguration) -> Material:
        # TODO: use create_slab_with_additional_layers, create_slab_with_additional_layers_float
        return create_slab_with_additional_layers(
            slab=configuration.crystal,
            number_of_additional_layers=configuration.number_of_additional_layers,
            vacuum=configuration.vacuum,
        )

    def _generate(self, configuration: SlabDefectConfiguration) -> Material:
        slab_with_additional_layers = self.get_slab_with_additional_layers(configuration)
        filtered_slab_with_additional_layers = filter_by_box(
            slab_with_additional_layers, max_coordinate=[1, 1, new_max_z], reset_ids=True
        )
