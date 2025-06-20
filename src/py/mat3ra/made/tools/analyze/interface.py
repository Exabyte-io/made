from functools import cached_property
from typing import List, Tuple, Dict

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.made.tools.build.slab.builders import SlabBuilder
from mat3ra.made.tools.build.slab.configuration import (
    SlabConfiguration,
    SlabStrainedSupercellConfiguration,
)
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from pymatgen.analysis.interfaces.coherent_interfaces import ZSLGenerator
from pymatgen.analysis.interfaces.zsl import vec_area


def compute_strains_from_match(
    substrate_sl: np.ndarray,
    film_sl: np.ndarray,
    substrate_orig: np.ndarray,
    film_orig: np.ndarray,
) -> Dict[str, float]:
    """
    Compute average, a-, and b-strains for substrate and film from superlattice vs original vectors.
    """

    def pair_strain(sl_vecs: np.ndarray, orig_vecs: np.ndarray) -> Tuple[float, float, float]:
        a = (np.linalg.norm(sl_vecs[0]) - np.linalg.norm(orig_vecs[0])) / np.linalg.norm(orig_vecs[0])
        b = (np.linalg.norm(sl_vecs[1]) - np.linalg.norm(orig_vecs[1])) / np.linalg.norm(orig_vecs[1])
        avg = (abs(a) + abs(b)) / 2
        return a, b, avg

    sa, sb, savg = pair_strain(substrate_sl, substrate_orig)
    fa, fb, favg = pair_strain(film_sl, film_orig)

    return {
        "substrate_strain": savg,
        "film_strain": favg,
        "substrate_a_strain": sa,
        "substrate_b_strain": sb,
        "film_a_strain": fa,
        "film_b_strain": fb,
    }


class BaseInterfaceAnalyzer(InMemoryEntityPydantic):
    substrate_slab_configuration: SlabConfiguration
    film_slab_configuration: SlabConfiguration

    @cached_property
    def substrate_material(self):
        return SlabBuilder().get_material(self.substrate_slab_configuration)

    @cached_property
    def film_material(self):
        return SlabBuilder().get_material(self.film_slab_configuration)

    def to_2d(self, lattice) -> np.ndarray:
        v = np.array(lattice.vector_arrays)
        return v[:2, :2]

    def to_3d_strain(self, strain2d: np.ndarray) -> Matrix3x3Schema:
        m = np.eye(3)
        m[:2, :2] = strain2d
        return Matrix3x3Schema(root=m.tolist())

    @property
    def identity_supercell(self) -> SupercellMatrix2DSchema:
        return SupercellMatrix2DSchema(root=[[1, 0], [0, 1]])

    @property
    def identity_strain(self) -> Matrix3x3Schema:
        return Matrix3x3Schema(root=np.eye(3).tolist())


class InterfaceAnalyzer(BaseInterfaceAnalyzer):
    """Simple interface analyzer without ZSL: computes direct strain and identity supercell."""

    @cached_property
    def film_strain_matrix(self) -> Matrix3x3Schema:
        sub2d = self.to_2d(self.substrate_material.lattice)
        film2d = self.to_2d(self.film_material.lattice)
        try:
            inv_film = np.linalg.inv(film2d)
        except np.linalg.LinAlgError:
            raise ValueError("Film lattice vectors are not linearly independent.")
        strain2d = inv_film @ sub2d
        return self.to_3d_strain(strain2d)

    @property
    def substrate_strain_matrix(self) -> Matrix3x3Schema:
        return self.identity_strain

    @property
    def substrate_supercell_matrix(self) -> SupercellMatrix2DSchema:
        return self.identity_supercell

    @property
    def film_supercell_matrix(self) -> SupercellMatrix2DSchema:
        return self.identity_supercell

    @cached_property
    def substrate_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return SlabStrainedSupercellConfiguration(
            stack_components=self.substrate_slab_configuration.stack_components,
            direction=self.substrate_slab_configuration.direction,
            xy_supercell_matrix=self.identity_supercell,
            strain_matrix=self.identity_strain,
        )

    @cached_property
    def film_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return SlabStrainedSupercellConfiguration(
            stack_components=self.film_slab_configuration.stack_components,
            direction=self.film_slab_configuration.direction,
            xy_supercell_matrix=self.identity_supercell,
            strain_matrix=self.film_strain_matrix,
        )


class ZSLConfigurationWithStrain(InMemoryEntityPydantic):
    """Configuration with ZSL strain metadata."""

    substrate_config: SlabStrainedSupercellConfiguration
    film_config: SlabStrainedSupercellConfiguration
    strain_info: Dict[str, float]


class ZSLInterfaceAnalyzer(BaseInterfaceAnalyzer):
    """Interface analyzer using Pymatgen's ZSL algorithm to find matching supercells."""

    max_area: float = 50.0
    max_area_ratio_tol: float = 0.09
    max_length_tol: float = 0.03
    max_angle_tol: float = 0.01

    @cached_property
    def _zsl(self) -> ZSLGenerator:
        return ZSLGenerator(
            max_area=self.max_area,
            max_area_ratio_tol=self.max_area_ratio_tol,
            max_length_tol=self.max_length_tol,
            max_angle_tol=self.max_angle_tol,
        )

    def get_strained_slab_configurations(
        self,
    ) -> List[Tuple[SlabStrainedSupercellConfiguration, SlabStrainedSupercellConfiguration]]:
        sub2d = self.to_2d(self.substrate_material.lattice)
        film2d = self.to_2d(self.film_material.lattice)
        sub_area = vec_area(sub2d[0], sub2d[1])
        film_area = vec_area(film2d[0], film2d[1])
        configs: List[Tuple[SlabStrainedSupercellConfiguration, SlabStrainedSupercellConfiguration]] = []

        sets = self._zsl.generate_sl_transformation_sets(film_area, sub_area)
        matches = self._zsl.get_equiv_transformations(sets, film2d, sub2d)

        for film_sl, sub_sl, film_trans, sub_trans in matches:
            # compute 2D strain
            sub_strain2d = sub_sl @ np.linalg.inv(sub2d @ sub_trans)
            film_strain2d = film_sl @ np.linalg.inv(film2d @ film_trans)

            sub_conf = SlabStrainedSupercellConfiguration(
                stack_components=self.substrate_slab_configuration.stack_components,
                direction=self.substrate_slab_configuration.direction,
                xy_supercell_matrix=SupercellMatrix2DSchema(root=sub_trans.astype(int).tolist()),
                strain_matrix=self.to_3d_strain(sub_strain2d),
            )
            film_conf = SlabStrainedSupercellConfiguration(
                stack_components=self.film_slab_configuration.stack_components,
                direction=self.film_slab_configuration.direction,
                xy_supercell_matrix=SupercellMatrix2DSchema(root=film_trans.astype(int).tolist()),
                strain_matrix=self.to_3d_strain(film_strain2d),
            )
            configs.append((sub_conf, film_conf))

        return configs

    def get_strained_slab_configurations_with_metadata(
        self,
    ) -> List[ZSLConfigurationWithStrain]:
        sub2d = self.to_2d(self.substrate_material.lattice)
        film2d = self.to_2d(self.film_material.lattice)
        sub_area = vec_area(sub2d[0], sub2d[1])
        film_area = vec_area(film2d[0], film2d[1])
        results: List[ZSLConfigurationWithStrain] = []

        sets = self._zsl.generate_sl_transformation_sets(film_area, sub_area)
        matches = self._zsl.get_equiv_transformations(sets, film2d, sub2d)

        for film_sl, sub_sl, film_trans, sub_trans in matches:
            sub_strain2d = sub_sl @ np.linalg.inv(sub2d @ sub_trans)
            film_strain2d = film_sl @ np.linalg.inv(film2d @ film_trans)

            sub_conf = SlabStrainedSupercellConfiguration(
                stack_components=self.substrate_slab_configuration.stack_components,
                direction=self.substrate_slab_configuration.direction,
                xy_supercell_matrix=SupercellMatrix2DSchema(root=sub_trans.astype(int).tolist()),
                strain_matrix=self.to_3d_strain(sub_strain2d),
            )
            film_conf = SlabStrainedSupercellConfiguration(
                stack_components=self.film_slab_configuration.stack_components,
                direction=self.film_slab_configuration.direction,
                xy_supercell_matrix=SupercellMatrix2DSchema(root=film_trans.astype(int).tolist()),
                strain_matrix=self.to_3d_strain(film_strain2d),
            )
            info = compute_strains_from_match(sub_sl, film_sl, sub2d, film2d)
            results.append(
                ZSLConfigurationWithStrain(
                    substrate_config=sub_conf,
                    film_config=film_conf,
                    strain_info=info,
                )
            )

        # sort by average total strain
        results.sort(key=lambda x: (x.strain_info["substrate_strain"] + x.strain_info["film_strain"]) / 2)
        return results
