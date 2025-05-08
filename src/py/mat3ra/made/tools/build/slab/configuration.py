from typing import List, Tuple

import numpy as np
from mat3ra.esse.models.materials_category.single_material.two_dimensional.slab.configuration import (
    SlabConfigurationSchema,
)


from mat3ra.made.material import Material

from .. import BaseConfiguration
from ...third_party import PymatgenSpacegroupAnalyzer
from ...convert import to_pymatgen, from_pymatgen


class SlabConfiguration(SlabConfigurationSchema, BaseConfiguration):
    """
    Configuration for building a slab.

    Args:
        bulk (Material): The bulk material.
        miller_indices (Tuple[int, int, int]): The Miller indices of the slab.
        thickness (int): The thickness of the slab.
        vacuum (float): The vacuum thickness, in Angstroms.
        xy_supercell_matrix (List[List[int]]): The supercell matrix for the xy plane.
        use_conventional_cell (bool): Whether to use the conventional cell.
        use_orthogonal_z (bool): Whether to use orthogonal z.
        make_primitive (bool): Whether to try to find primitive cell for the created slab.
    """

    bulk: Material

    def __init__(
        self,
        bulk=None,
        miller_indices=SlabConfigurationSchema.model_fields["miller_indices"].default,
        thickness=SlabConfigurationSchema.model_fields["thickness"].default,
        vacuum=SlabConfigurationSchema.model_fields["vacuum"].default,
        xy_supercell_matrix=None,
        use_conventional_cell=SlabConfigurationSchema.model_fields["use_conventional_cell"].default,
        use_orthogonal_z=SlabConfigurationSchema.model_fields["use_orthogonal_z"].default,
        make_primitive=SlabConfigurationSchema.model_fields["make_primitive"].default,
    ):
        if xy_supercell_matrix is None:
            xy_supercell_matrix = np.eye(2).tolist()
        bulk = bulk or Material.create_default()
        __bulk_pymatgen_structure = (
            PymatgenSpacegroupAnalyzer(to_pymatgen(bulk)).get_conventional_standard_structure()
            if use_conventional_cell
            else to_pymatgen(bulk)
        )
        __bulk_config = from_pymatgen(__bulk_pymatgen_structure)
        super().__init__(
            bulk=Material.create(__bulk_config),
            miller_indices=miller_indices,
            thickness=thickness,
            vacuum=vacuum,
            xy_supercell_matrix=xy_supercell_matrix,
            use_conventional_cell=use_conventional_cell,
            use_orthogonal_z=use_orthogonal_z,
            make_primitive=make_primitive,
        )

    @property
    def _json(self):
        return {
            "type": "SlabConfiguration",
            "bulk": self.bulk.to_dict(),
            "miller_indices": self.miller_indices,
            "thickness": self.thickness,
            "vacuum": self.vacuum,
            "xy_supercell_matrix": self.xy_supercell_matrix,
            "use_conventional_cell": self.use_conventional_cell,
            "use_orthogonal_z": self.use_orthogonal_z,
            "make_primitive": self.make_primitive,
        }
