from typing import List, Union

from mat3ra.made.material import Material
from ..defect.point.builders import AtomAtCoordinateConfiguration
from ..merge.configuration import MergeConfiguration


class PassivationConfiguration(MergeConfiguration):
    """
    Configuration for a passivation.

    Args:
        merge_components (List[Material, AtomAtCoordinateConfiguration]): The components to merge.
        passivant (str): The passivating element.
        bond_length (float): The bond length.
    """

    merge_components: List[Union[Material, AtomAtCoordinateConfiguration]]
    passivant: str = "H"
    bond_length: float = 1.0
