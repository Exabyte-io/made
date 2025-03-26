import numpy as np
from mat3ra.esse.models.materials_category.single_material.two_dimensional.slab.configuration import (
    SlabConfigurationSchema,
)

from mat3ra.made.material import Material

from .. import BaseConfiguration
from ...third_party import PymatgenSpacegroupAnalyzer
from ...convert import to_pymatgen, from_pymatgen


class SlabConfiguration(BaseConfiguration, SlabConfigurationSchema):
    def __init__(
        self,
        bulk=SlabConfigurationSchema.bulk,
        miller_indices=SlabConfigurationSchema.miller_indices,
        thickness=SlabConfigurationSchema.thickness,
        vacuum=SlabConfigurationSchema.vacuum,
        xy_supercell_matrix=SlabConfigurationSchema.xy_supercell_matrix,
        use_conventional_cell=SlabConfigurationSchema.use_conventional_cell,
        use_orthogonal_z=SlabConfigurationSchema.use_orthogonal_z,
        make_primitive=SlabConfigurationSchema.make_primitive,
    ):
        if xy_supercell_matrix is None:
            xy_supercell_matrix = np.eye(2).tolist()
        bulk = bulk or Material(Material.default_config)
        __bulk_pymatgen_structure = (
            PymatgenSpacegroupAnalyzer(to_pymatgen(bulk)).get_conventional_standard_structure()
            if use_conventional_cell
            else to_pymatgen(bulk)
        )
        __bulk_config = from_pymatgen(__bulk_pymatgen_structure)
        super().__init__(
            bulk=Material(__bulk_config),
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
            "bulk": self.bulk.to_json(),
            "miller_indices": self.miller_indices,
            "thickness": self.thickness,
            "vacuum": self.vacuum,
            "xy_supercell_matrix": self.xy_supercell_matrix,
            "use_conventional_cell": self.use_conventional_cell,
            "use_orthogonal_z": self.use_orthogonal_z,
            "make_primitive": self.make_primitive,
        }
