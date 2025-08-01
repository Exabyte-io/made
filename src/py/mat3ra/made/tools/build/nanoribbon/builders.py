from typing import Any, Optional, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.modify import translate_to_center
from . import NanoribbonConfiguration
from .. import MaterialWithBuildMetadata
from ..nanotape import NanoTapeConfiguration
from ..nanotape.builders import NanoTapeBuilder, NanoTapeBuilderParameters


class NanoribbonBuilderParameters(NanoTapeBuilderParameters):
    pass


class NanoribbonBuilder(NanoTapeBuilder):
    _ConfigurationType = "NanoribbonConfiguration"  # String type annotation to avoid circular import
    _BuilderParametersType = NanoribbonBuilderParameters
    _DefaultBuildParameters = NanoribbonBuilderParameters(
        use_rectangular_lattice=True,
    )

    @property
    def stack_component_types_conversion_map(self):
        return {**super().stack_component_types_conversion_map, NanoTapeConfiguration: NanoTapeBuilder}

    def _post_process(self, item: Material, post_process_parameters: Optional[Any] = None) -> Material:
        item = super()._post_process(item, post_process_parameters)
        item = translate_to_center(item, axes=["x", "y"])
        return item

    def _update_material_name(
        self, material: Union[Material, MaterialWithBuildMetadata], configuration: Any
    ) -> Material:
        if isinstance(configuration, NanoribbonConfiguration):
            nanotape = configuration.nanotape
            material = self._update_material_name_with_edge_type(
                material, nanotape.lattice_lines.crystal.name, nanotape.lattice_lines.miller_indices_2d, "Nanoribbon"
            )
            return material
        return super()._update_material_name(material, configuration)

    def _update_material_name_with_edge_type(
        self,
        material: Union[Material, MaterialWithBuildMetadata],
        crystal_name: str,
        miller_indices_2d: tuple,
        structure_type: str,
    ) -> Material:
        edge_type = self._get_edge_type_from_miller_indices(miller_indices_2d)
        miller_str = f"{miller_indices_2d[0]}{miller_indices_2d[1]}"
        material.name = f"{crystal_name} - {edge_type} {structure_type} ({miller_str})"
        return material
