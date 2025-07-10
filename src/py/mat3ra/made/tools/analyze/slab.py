from mat3ra.code.entity import InMemoryEntityPydantic

from ..analyze import BaseMaterialAnalyzer
from ..build.metadata import MaterialMetadata
from ..build.slab.configurations import SlabConfiguration
from ..build.vacuum.configuration import VacuumConfiguration


class SlabMaterialAnalyzer(BaseMaterialAnalyzer):
    def get_slab_configuration(self) -> SlabConfiguration:
        metadata = MaterialMetadata(**self.material.metadata)
        slab_build_configuration_dict = metadata.build[-1].configuration
        if slab_build_configuration_dict.get("type") != "SlabConfiguration":
            raise ValueError("Material is not a slab.")

        return SlabConfiguration(**slab_build_configuration_dict)

    @property
    def number_of_layers(self) -> int:
        slab_configuration = self.get_slab_configuration()
        return slab_configuration.number_of_layers

    @property
    def vacuum_ratio(self) -> float:
        slab_configuration = self.get_slab_configuration()
        return slab_configuration.vacuum / self.material.lattice.c

    @property
    def vacuum_thickness_in_layers(self) -> float:
        slab_configuration = self.get_slab_configuration()
        return self.vacuum_ratio / (1 - self.vacuum_ratio) * self.number_of_layers

    def get_slab_configuration_with_no_vacuum(self) -> SlabConfiguration:
        """
        Get the slab material without vacuum.
        Returns:
            InMemoryEntityPydantic: The slab material without vacuum.
        """
        slab_configuration = self.get_slab_configuration()
        slab_configuration_with_no_vacuum = slab_configuration.clone()
        slab_configuration_with_no_vacuum.set_vacuum(0.0)

        return slab_configuration_with_no_vacuum

    def get_slab_vacuum_configuration(self) -> VacuumConfiguration:
        """
        Get the vacuum configuration from the slab.
        Returns:
            VacuumConfiguration: The vacuum configuration.
        """
        slab_configuration = self.get_slab_configuration()
        return slab_configuration.vacuum_configuration
