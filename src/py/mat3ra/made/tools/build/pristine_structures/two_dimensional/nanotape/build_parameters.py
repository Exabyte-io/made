from pydantic import Field

from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseBuilderParameters


class NanoTapeBuilderParameters(BaseBuilderParameters):
    use_rectangular_lattice: bool = Field(
        True, description="If True, set the XY lattice to be rectangular after stacking."
    )
