from ...convert import from_ase
from ...third_party import ASEAtoms


class ConvertGeneratedItemsASEAtomsMixin:
    @staticmethod
    def _convert_generated_item(item: ASEAtoms):
        return from_ase(item)
