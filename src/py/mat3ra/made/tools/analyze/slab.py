from ..analyze import BaseMaterialAnalyzer
from ..build.slab.configurations import SlabConfiguration
from ..build.slab.build_parameters import SlabBuilderParameters
from ..build.vacuum.configuration import VacuumConfiguration


class SlabMaterialAnalyzer(BaseMaterialAnalyzer):
    def _traverse_build_history_for_slab(self):
        """
        Traverse build history from the end to find slab configuration and parameters.
        
        Returns:
            tuple: (slab_configuration_dict, build_parameters_dict, supercell_matrix)
        """
        supercell_matrix = None
        build_parameters_dict = None
        slab_configuration_dict = None

        for build_step in reversed(self.material.metadata["build"]):
            config = build_step["configuration"]
            config_type = config.get("type")
            if config_type == "SlabConfiguration":
                slab_configuration_dict = config
                build_parameters_dict = build_step.get("build_parameters")
                break
            elif config_type == "SupercellConfiguration":
                supercell_matrix = config.get("supercell_matrix")

        if slab_configuration_dict is None:
            raise ValueError("Material is not a slab.")

        return slab_configuration_dict, build_parameters_dict, supercell_matrix

    def get_slab_configuration(self) -> SlabConfiguration:
        slab_configuration_dict, _, _ = self._traverse_build_history_for_slab()
        return SlabConfiguration(**slab_configuration_dict)

    def get_build_parameters(self) -> SlabBuilderParameters:
        _, build_parameters_dict, supercell_matrix = self._traverse_build_history_for_slab()

        if supercell_matrix:
            if build_parameters_dict:
                build_parameters_dict["xy_supercell_matrix"] = supercell_matrix
            else:
                build_parameters_dict = {"xy_supercell_matrix": supercell_matrix}

        return SlabBuilderParameters(**build_parameters_dict) if build_parameters_dict else SlabBuilderParameters()

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
        return self.vacuum_ratio / (1 - self.vacuum_ratio) * self.number_of_layers

    def get_slab_configuration_with_no_vacuum(self) -> SlabConfiguration:
        slab_configuration = self.get_slab_configuration()
        slab_configuration_with_no_vacuum = slab_configuration.clone()
        slab_configuration_with_no_vacuum.set_vacuum(0.0)

        return slab_configuration_with_no_vacuum

    def get_slab_vacuum_configuration(self) -> VacuumConfiguration:
        slab_configuration = self.get_slab_configuration()
        return slab_configuration.vacuum_configuration
