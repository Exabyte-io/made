from mat3ra.esse.models.materials_category.single_material.two_dimensional.slab.configuration import (
    SlabConfigurationSchema,
)
from mat3ra.made.material import Material
from .. import BaseConfiguration
from ...third_party import PymatgenSpacegroupAnalyzer
from ...convert import to_pymatgen, from_pymatgen


class SlabConfiguration(SlabConfigurationSchema, BaseConfiguration):
    def __init__(self, **data):
        super().__init__(**data)
        if self.bulk is None:
            self.bulk = Material(Material.default_config)

        bulk_structure = to_pymatgen(self.bulk)
        if self.use_conventional_cell:
            bulk_structure = PymatgenSpacegroupAnalyzer(bulk_structure).get_conventional_standard_structure()

        self.bulk = Material(from_pymatgen(bulk_structure))

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
