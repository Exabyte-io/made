from typing import Optional, Any, Union, Tuple

from mat3ra.made.material import Material
from .....build_components import MaterialWithBuildMetadata, TypeConfiguration
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.builder import SlabBuilder
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.configuration import SlabConfiguration
from .....modify import translate_to_z_level, filter_by_box
from .configurations import MonolayerConfiguration


class MonolayerBuilder(SlabBuilder):
    """
    Builder for creating monolayer structures from crystal materials.

    The builder creates different monolayer structures based on the crystal type:
    - HEX: Creates a slab with Miller indices (0,0,1), thickness=1, then filters half
    - FCC/CUB: Creates a slab with Miller indices (1,1,1), thickness=1 using primitive cell
    """

    _ConfigurationType = MonolayerConfiguration
    _GeneratedItemType = Material

    def _get_miller_indices_for_lattice_type(self, crystal) -> Tuple[int, int, int]:
        lattice_type_str = str(crystal.lattice.type.value)

        if lattice_type_str == "HEX":
            return (0, 0, 1)
        else:
            return (1, 1, 1)

    def _generate(self, configuration: MonolayerConfiguration) -> Material:
        crystal = configuration.crystal
        miller_indices = self._get_miller_indices_for_lattice_type(crystal)
        vacuum = configuration.vacuum

        slab_config = SlabConfiguration.from_parameters(
            material_or_dict=crystal,
            miller_indices=miller_indices,
            number_of_layers=1,
            vacuum=vacuum,
            use_conventional_cell=False,
        )
        slab = super()._generate(slab_config)

        lattice_type_str = crystal.lattice.type.value
        if lattice_type_str == "HEX":
            slab = self._apply_hex_filtering(slab)

        return slab

    def _apply_hex_filtering(self, slab: Material) -> Material:
        centered_slab = translate_to_z_level(slab, z_level="center")

        half_filtered_slab = filter_by_box(
            centered_slab, min_coordinate=[0, 0, 0.0], max_coordinate=[1, 1, 0.5], use_cartesian_coordinates=False
        )

        return half_filtered_slab

    # TODO: we need to move this to a common place for this and CLPBuilder
    def _enforce_convention(self, material: Union[Material, MaterialWithBuildMetadata]) -> Material:
        return translate_to_z_level(material, "bottom")

    def _post_process(
        self,
        item: Material,
        post_process_parameters: Optional[Any] = None,
        configuration: Optional[TypeConfiguration] = None,
    ) -> Material:
        item = super()._post_process(item, post_process_parameters, configuration)
        return self._enforce_convention(item)

    def _update_material_name(
        self, material: Union[Material, MaterialWithBuildMetadata], configuration: MonolayerConfiguration
    ) -> Material:
        crystal = configuration.crystal

        original_name = crystal.name or "Crystal"
        material.name = f"{original_name} - Monolayer"
        return material
