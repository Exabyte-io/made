import numpy as np
from typing import Union, List, Optional, Tuple
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab_strained_supercell import (
    SlabStrainedSupercellConfigurationSchema,
)

from mat3ra.made.material import Material
from mat3ra.code.constants import AtomicCoordinateUnits
from .slab_configuration import SlabConfiguration


class SlabStrainedSupercellConfiguration(SlabConfiguration, SlabStrainedSupercellConfigurationSchema):
    type: str = "SlabStrainedSupercellConfiguration"
    strain_matrix: Matrix3x3Schema = Matrix3x3Schema(root=np.eye(3).tolist())

    @classmethod
    def from_parameters(
        cls,
        material_or_dict: Union[Material, dict],
        miller_indices: Tuple[int, int, int],
        number_of_layers: int,
        termination_formula: Optional[str] = None,
        vacuum: float = 10.0,
        use_conventional_cell: bool = True,
        strain_matrix: Optional[List[List[float]]] = None,
    ) -> "SlabStrainedSupercellConfiguration":
        config = super().from_parameters(
            material_or_dict=material_or_dict,
            miller_indices=miller_indices,
            number_of_layers=number_of_layers,
            termination_formula=termination_formula,
            vacuum=vacuum,
            use_conventional_cell=use_conventional_cell,
        )
        if strain_matrix:
            config.strain_matrix = Matrix3x3Schema(root=strain_matrix)
        if isinstance(material_or_dict, dict):
            config.atomic_layers.crystal.basis.units = AtomicCoordinateUnits.crystal
        return config
