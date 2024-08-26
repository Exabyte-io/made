from mat3ra.made.material import Material

from .enums import SurfaceTypes, EdgeTypes
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

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "slab": self.slab.to_json(),
            "passivant": self.passivant,
            "bond_length": self.bond_length,
        }


class SurfacePassivationConfiguration(PassivationConfiguration):
    """
    Configuration for a passivation.
    """

    surface: SurfaceTypes = SurfaceTypes.BOTH


class EdgePassivationConfiguration(PassivationConfiguration):
    """
    Configuration for a passivation.
    """

    edge: EdgeTypes = EdgeTypes.BOTH
