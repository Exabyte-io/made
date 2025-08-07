from .configuration import AdatomDefectConfiguration
from .....build_components.entities.reusable.two_dimensional.slab_stack.builder import SlabStackBuilder


class AdatomDefectBuilder(SlabStackBuilder):
    _ConfigurationType = AdatomDefectConfiguration

    def get_name_suffix(self, configuration: AdatomDefectConfiguration) -> str:
        elements = set(configuration.added_component.basis.elements.values)
        elements_str = ", ".join(sorted(elements))
        return f"{elements_str} Adatom"
