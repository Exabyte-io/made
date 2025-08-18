from ...convert import from_pymatgen
from ...third_party import PymatgenStructure


class ConvertGeneratedItemsPymatgenStructureMixin:
    @staticmethod
    def _convert_generated_item(item: PymatgenStructure):
        return from_pymatgen(item)
