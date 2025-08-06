from typing import Union, List, Optional

import numpy as np
from mat3ra.esse.models.materials_category.compound_pristine_structures.two_dimensional.interface.configuration import (  # noqa: E501
    InterfaceConfigurationSchema,
)

from .....pristine_structures.two_dimensional.slab_strained_supercell.configuration import (
    SlabStrainedSupercellConfiguration,
)
from ......analyze.utils import calculate_von_mises_strain
from ......build_components.entities.core.two_dimensional.vacuum.configuration import VacuumConfiguration
from ......build_components.operations.core.combinations.stack.configuration import StackConfiguration
from ......utils import unwrap


class InterfaceConfiguration(StackConfiguration, InterfaceConfigurationSchema):
    # components and their modifiers added in the order they are stacked, from bottom to top
    stack_components: List[
        Union[
            SlabStrainedSupercellConfiguration,
            VacuumConfiguration,
        ]
    ]
    xy_shift: List[float] = InterfaceConfigurationSchema.model_fields["xy_shift"].default  # in Angstroms

    type: str = "InterfaceConfiguration"

    @property
    def substrate_configuration(self) -> SlabStrainedSupercellConfiguration:
        return self.stack_components[0]

    @property
    def film_configuration(self) -> SlabStrainedSupercellConfiguration:
        return self.stack_components[1]

    @property
    def vacuum_configuration(self) -> Optional[VacuumConfiguration]:
        if len(self.stack_components) > 2:
            return self.stack_components[2]
        return None

    @property
    def von_mises_strain_percentage(self) -> float:
        """Von Mises strain (%) from the film's strain_matrix (2D part)."""
        raw = self.film_configuration.strain_matrix.root
        strain_matrix = np.array(unwrap(raw))
        return calculate_von_mises_strain(strain_matrix)
