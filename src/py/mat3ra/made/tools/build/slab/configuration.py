from mat3ra.esse.models.materials_category.single_material.two_dimensional.slab.configuration import (
    SlabConfigurationSchema,
)

from mat3ra.made.material import Material
from .. import BaseConfigurationPydantic
from ...convert import to_pymatgen, from_pymatgen
from ...third_party import PymatgenSpacegroupAnalyzer


class SlabConfiguration(SlabConfigurationSchema, BaseConfigurationPydantic):
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

    type: str = "SlabConfiguration"
    bulk: Material

    def __init__(self, **kwargs):
        bulk = kwargs.pop("bulk", None) or Material.create_default()
        use_conventional = kwargs.get(
            "use_conventional_cell", SlabConfigurationSchema.model_fields["use_conventional_cell"].default
        )

        bulk_pmg = to_pymatgen(bulk)
        if use_conventional:
            bulk_pmg = PymatgenSpacegroupAnalyzer(bulk_pmg).get_conventional_standard_structure()

        bulk_config = from_pymatgen(bulk_pmg)
        kwargs["bulk"] = Material.create(bulk_config)

        super().__init__(**kwargs)
