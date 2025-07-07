from pydantic import BaseModel

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.slab.configuration import SlabDefectConfiguration
from mat3ra.made.tools.build.merge import MergeBuilder


class SlabDefectBuilderParameters(BaseModel):
    auto_add_vacuum: bool = True
    vacuum_thickness: float = 5.0


class SlabDefectBuilder(MergeBuilder):
    _ConfigurationType = SlabDefectConfiguration

    def _generate(self, configuration: SlabDefectConfiguration) -> Material:
        return super()._generate(configuration)
