from typing import Any

from mat3ra.made.material import Material
from . import NanoribbonConfiguration
from ..nanotape.builders import NanoTapeBuilder, NanoTapeBuilderParameters


class NanoribbonBuilderParameters(NanoTapeBuilderParameters):
    pass


class NanoribbonBuilder(NanoTapeBuilder):
    _ConfigurationType = "NanoribbonConfiguration"  # String type annotation to avoid circular import
    _BuilderParametersType = NanoribbonBuilderParameters
    _DefaultBuildParameters = NanoribbonBuilderParameters(
        use_rectangular_lattice=True,
    )

    def _update_material_name(self, material: Material, configuration: Any) -> Material:
        if isinstance(configuration, NanoribbonConfiguration):
            nanotape = configuration.nanotape
            material = self._update_material_name_with_edge_type(
                material, nanotape.lattice_lines.crystal.name, nanotape.lattice_lines.miller_indices_2d, "Nanoribbon"
            )
            return material
        return super()._update_material_name(material, configuration)

    def _update_material_name_with_edge_type(
        self, material: Material, crystal_name: str, miller_indices_2d: tuple, structure_type: str
    ) -> Material:
        edge_type = self._get_edge_type_from_miller_indices(miller_indices_2d)
        miller_str = f"{miller_indices_2d[0]}{miller_indices_2d[1]}"
        material.name = f"{crystal_name} - {edge_type} {structure_type} ({miller_str})"
        return material
