from mat3ra.made.material import Material

from ...enums import SurfaceTypes
from ...build import BaseConfiguration


class PassivationConfiguration(BaseConfiguration):
    """
    Configuration for a passivation.

    Args:
        slab (Material): The Material object.
        passivant (str): The passivating element.
        bond_length (float): The bond length.
    """

    slab: Material
    passivant: str = "H"
    bond_length: float = 1.0
    surface: SurfaceTypes = SurfaceTypes.BOTH

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "slab": self.slab.to_json(),
            "passivant": self.passivant,
            "bond_length": self.bond_length,
            "surface": self.surface.value,
        }
