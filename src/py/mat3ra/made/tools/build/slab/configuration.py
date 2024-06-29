from typing import List, Tuple, Any
from pydantic import BaseModel

from mat3ra.code.entity import InMemoryEntity

from mat3ra.made.material import Material
from ...third_party import PymatgenSpacegroupAnalyzer
from ...convert import to_pymatgen, from_pymatgen


class SlabConfiguration(BaseModel, InMemoryEntity):
    """
    Configuration for building a slab.

    Args:
        bulk (Material): The bulk material.
        miller_indices (Tuple[int, int, int]): The Miller indices of the slab.
        thickness (int): The thickness of the slab.
        vacuum (float): The vacuum thickness.
        xy_supercell_matrix (List[List[int]]): The supercell matrix for the xy plane.
        use_conventional_cell (bool): Whether to use the conventional cell.
        use_orthogonal_z (bool): Whether to use orthogonal z.
    """

    # TODO: fix arbitrary_types_allowed error and set Material class type
    bulk: Any
    miller_indices: Tuple[int, int, int] = (0, 0, 1)
    thickness: int = 1
    vacuum: int = 1
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]
    use_conventional_cell: bool = True
    use_orthogonal_z: bool = False

    def __init__(
        self,
        bulk=None,
        miller_indices=miller_indices,
        thickness=thickness,
        vacuum=vacuum,
        xy_supercell_matrix=xy_supercell_matrix,
        use_conventional_cell=use_conventional_cell,
        use_orthogonal_z=use_orthogonal_z,
    ):
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
            use_orthogonal_z=use_orthogonal_z,
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
        }
