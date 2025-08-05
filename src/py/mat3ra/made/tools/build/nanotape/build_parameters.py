from pydantic import Field

from .. import BaseBuilderParameters


class NanoTapeBuilderParameters(BaseBuilderParameters):
    use_rectangular_lattice: bool = Field(
        True, description="If True, set the XY lattice to be rectangular after stacking."
    )
