from mat3ra.made.material import Material
from .configurations import MonolayerConfiguration
from ..slab.builders import SlabBuilder
from ..slab.configurations import SlabConfiguration
from ...modify import translate_to_z_level, filter_by_box


class MonolayerBuilder(SlabBuilder):
    """
    Builder for creating monolayer structures from crystal materials.

    The builder creates different monolayer structures based on the crystal type:
    - HEX: Creates a slab with Miller indices (0,0,1), thickness=1, then filters half
    - FCC/CUB: Creates a slab with Miller indices (1,1,1), thickness=1 using primitive cell
    """

    _ConfigurationType = MonolayerConfiguration
    _GeneratedItemType = Material

    def _generate(self, configuration: MonolayerConfiguration) -> Material:
        crystal = configuration.crystal
        lattice_type_str = (
            crystal.lattice.type.value if hasattr(crystal.lattice.type, "value") else str(crystal.lattice.type)
        )
        vacuum = configuration.vacuum

        if lattice_type_str == "HEX":
            miller_indices = (0, 0, 1)
        elif lattice_type_str in ["FCC", "CUB"]:
            miller_indices = (1, 1, 1)
        else:
            miller_indices = (1, 1, 1)

        slab_config = SlabConfiguration.from_parameters(
            material_or_dict=crystal,
            miller_indices=miller_indices,
            number_of_layers=1,
            vacuum=vacuum,
        )
        slab = super()._generate(slab_config)

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
    def _enforce_convention(self, material: Material) -> Material:
        return translate_to_z_level(material, "bottom")

    def _post_process(self, item: Material, post_process_parameters=None) -> Material:
        item = super()._post_process(item, post_process_parameters)
        return self._enforce_convention(item)

    def _update_material_name(self, material: Material, configuration: MonolayerConfiguration) -> Material:
        crystal = configuration.crystal
        lattice_type_str = (
            crystal.lattice.type.value if hasattr(crystal.lattice.type, "value") else str(crystal.lattice.type)
        )

        if lattice_type_str == "HEX":
            miller_str = "001"
        elif lattice_type_str in ["FCC", "CUB"]:
            miller_str = "111"
        else:
            miller_str = "111"

        original_name = crystal.name or "Crystal"
        material.name = f"{original_name} - Monolayer ({miller_str})"
        return material
