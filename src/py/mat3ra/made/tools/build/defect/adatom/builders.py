from .configuration import AdatomDefectConfiguration
from ..slab.builders import SlabStackBuilder


class AdatomDefectBuilder(SlabStackBuilder):
    _ConfigurationType = AdatomDefectConfiguration

    def get_name_suffix(self, configuration: AdatomDefectConfiguration) -> str:
        elements = set(configuration.added_component.basis.elements.values)
        elements_str = ", ".join(sorted(elements))
        return f"{elements_str} Adatom"
