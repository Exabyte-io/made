from typing import List, Union

from mat3ra.esse.models.materials_category.processed_structures.two_dimensional.passivation.configuration import (
    PassivationConfigurationSchema,
)
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.material import Material
from .....build_components import MaterialWithBuildMetadata
from .....build_components.entities.core.zero_dimensional.atom.configuration import AtomAtCoordinateConfiguration
from .....build_components.operations.core.combinations.merge.configuration import MergeConfiguration


class PassivationConfiguration(MergeConfiguration, PassivationConfigurationSchema):
    """
    Configuration for a passivation.

    Args:
        merge_components (List[Material, AtomAtCoordinateConfiguration]): The components to merge.
        passivant (str): The passivating element.
        bond_length (float): The bond length.
    """

    merge_components: List[Union[Material, MaterialWithBuildMetadata, AtomAtCoordinateConfiguration]]
    merge_method: MergeMethodsEnum = MergeMethodsEnum.ADD
    passivant: str = "H"
    bond_length: float = 1.0

    @property
    def material(self) -> MaterialWithBuildMetadata:
        return self.merge_components[0]
