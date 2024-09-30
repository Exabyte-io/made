from typing import List, Tuple

import numpy as np

from mat3ra.made.material import Material

from .. import BaseConfiguration
from ...third_party import PymatgenSpacegroupAnalyzer
from ...convert import to_pymatgen, from_pymatgen


class SlabConfiguration(BaseConfiguration):
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
    miller_indices: Tuple[int, int, int] = (0, 0, 1)
    thickness: int = 1
    vacuum: float = 5.0
    xy_supercell_matrix: List[List[int]] = np.eye(2).tolist()
    use_conventional_cell: bool = True
    use_orthogonal_z: bool = False
    make_primitive: bool = False

    def __init__(
        self,
        bulk=None,
        miller_indices=miller_indices,
        thickness=thickness,
        vacuum=vacuum,
        xy_supercell_matrix=None,
        use_conventional_cell=use_conventional_cell,
        use_orthogonal_z=use_orthogonal_z,
        make_primitive=make_primitive,
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
