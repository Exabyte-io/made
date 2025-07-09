from pydantic import BaseModel

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.slab.configuration import SlabDefectConfiguration, AdatomDefectConfiguration
from mat3ra.made.tools.build.defect.enums import AdatomPlacementMethodEnum
from mat3ra.made.tools.build.merge import MergeBuilder
from mat3ra.made.tools.build.defect.point.builders import PointDefectBuilder


class SlabDefectBuilderParameters(BaseModel):
    auto_add_vacuum: bool = True
    vacuum_thickness: float = 5.0


class SlabDefectBuilder(MergeBuilder):
    _ConfigurationType = SlabDefectConfiguration

    def _generate(self, configuration: SlabDefectConfiguration) -> Material:
        return super()._generate(configuration)


class AdatomDefectBuilder(PointDefectBuilder):
    _ConfigurationType = AdatomDefectConfiguration
