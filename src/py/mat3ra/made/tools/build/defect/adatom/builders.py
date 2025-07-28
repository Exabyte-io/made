from mat3ra.made.tools.analyze import BaseMaterialAnalyzer
from .configuration import AdatomDefectConfiguration
from ..slab.builders import SlabStackBuilder


class AdatomDefectBuilder(SlabStackBuilder):
    _ConfigurationType = AdatomDefectConfiguration

    def get_name_suffix(self, configuration: AdatomDefectConfiguration) -> str:
        formula = BaseMaterialAnalyzer(material=configuration.added_component).formula
        return f"{formula} Adatom"
