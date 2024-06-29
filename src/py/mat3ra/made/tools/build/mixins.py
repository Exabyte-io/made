from ..convert import from_ase, from_pymatgen
from ..third_party import ASEAtoms, PymatgenStructure


class ConvertGeneratedItemsASEAtomsMixin:
    @staticmethod
    def _convert_generated_item(item: ASEAtoms):
        return from_ase(item)


class ConvertGeneratedItemsPymatgenStructureMixin:
    @staticmethod
    def _convert_generated_item(item: PymatgenStructure):
        return from_pymatgen(item)
