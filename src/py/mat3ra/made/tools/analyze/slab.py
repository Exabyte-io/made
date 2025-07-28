from typing import Type

from mat3ra.made.utils import adjust_material_cell_to_set_gap_along_direction, get_atomic_coordinates_extremum

from ..build import BuildMetadata, MaterialWithBuildMetadata
from ..build.slab.builders import SlabBuilderParameters
from ..build.slab.configurations import SlabConfiguration
from ..build.vacuum.configuration import VacuumConfiguration
from . import BaseMaterialAnalyzer
from .crystal_site import CrystalSiteAnalyzer


class BuildMetadataAnalyzer(BaseMaterialAnalyzer):
    material: MaterialWithBuildMetadata
    configuration_cls: None
    build_parameters_cls: None

    @property
    def build_metadata(self) -> BuildMetadata:
        return self.material.metadata.get_build_metadata_of_type(self.configuration_cls.__name__)

    @property
    def build_configuration(self) -> "configuration_cls":
        return self.configuration_cls(**self.build_metadata.configuration)

    @property
    def build_parameters(self) -> "configuration_cls":
        return self.build_parameters_cls(**self.build_metadata.build_parameters)


class SlabMaterialAnalyzer(BuildMetadataAnalyzer, CrystalSiteAnalyzer):
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
