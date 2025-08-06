from __future__ import annotations

from typing import Type

from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.build_parameters import SlabBuilderParameters
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.configuration import SlabConfiguration
from mat3ra.made.utils import adjust_material_cell_to_set_gap_along_direction, get_atomic_coordinates_extremum

from ..build_components.entities.core.two_dimensional.vacuum.configuration import VacuumConfiguration
from ..build_components.metadata import MaterialWithBuildMetadata
from .build_metadata_analyzer import BuildMetadataAnalyzer
from .crystal_site.crystal_site_analyzer import CrystalSiteAnalyzer


class SlabMaterialAnalyzer(
    BuildMetadataAnalyzer[SlabConfiguration, SlabBuilderParameters],
    CrystalSiteAnalyzer,
):
    configuration_cls: Type[SlabConfiguration] = SlabConfiguration
    build_parameters_cls: Type[SlabBuilderParameters] = SlabBuilderParameters

    @property
    def number_of_layers(self) -> int:
        return self.build_configuration.number_of_layers

    @property
    def layer_thickness(self) -> float:
        return (self.material.lattice.c - self.build_configuration.vacuum) / self.number_of_layers

    @property
    def vacuum_ratio(self) -> float:
        return self.build_configuration.vacuum / self.material.lattice.c

    @property
    def vacuum_thickness_in_layers(self) -> float:
        return self.vacuum_ratio / (1 - self.vacuum_ratio) * self.number_of_layers

    @property
    def slab_configuration_with_no_vacuum(self) -> SlabConfiguration:
        slab_configuration_with_no_vacuum = self.build_configuration.clone()
        slab_configuration_with_no_vacuum.set_vacuum(0.0)
        return slab_configuration_with_no_vacuum

    def get_slab_vacuum_configuration(self) -> VacuumConfiguration:
        return self.build_configuration.vacuum_configuration

    @property
    def max_z_crystal(self) -> float:
        return get_atomic_coordinates_extremum(self.material, "max", "z")

    @property
    def min_z_crystal(self) -> float:
        return get_atomic_coordinates_extremum(self.material, "min", "z")

    def get_gap_from_top(self, use_cartesian=False) -> float:
        gap_top_in_crystal = 1.0 - self.max_z_crystal
        if use_cartesian:
            gap_top_in_cartesian = self.material.basis.cell.convert_point_to_cartesian([0, 0, gap_top_in_crystal])[2]
            return gap_top_in_cartesian
        return gap_top_in_crystal

    def get_gap_from_bottom(self, use_cartesian=False) -> float:
        gap_bottom_in_crystal = self.min_z_crystal
        if use_cartesian:
            gap_bottom_in_cartesian = self.material.basis.cell.convert_point_to_cartesian(
                [0, 0, gap_bottom_in_crystal]
            )[2]
            return gap_bottom_in_cartesian
        return gap_bottom_in_crystal

    @property
    def _slab_with_no_gap(self) -> MaterialWithBuildMetadata:
        """
        Returns a fictitious structure with no gap along the z-direction. Used for stacking components.
        """
        return adjust_material_cell_to_set_gap_along_direction(self.material, 0)
