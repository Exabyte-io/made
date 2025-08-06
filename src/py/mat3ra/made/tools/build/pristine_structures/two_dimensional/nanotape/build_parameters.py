from pydantic import Field

from mat3ra.made.tools.build_components.entities.reusable.three_dimensional.crystal_lattice_base.build_parameters import \
    BaseBuilderParameters


class NanoTapeBuilderParameters(BaseBuilderParameters):
    use_rectangular_lattice: bool = Field(
        True, description="If True, set the XY lattice to be rectangular after stacking."
    )
