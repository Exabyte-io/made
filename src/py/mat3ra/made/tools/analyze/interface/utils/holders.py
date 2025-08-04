from typing import Optional

from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.made.tools.build.slab.strained_supercell_slab.configuration import SlabStrainedSupercellConfiguration


class MatchedSubstrateFilmConfigurationHolder(InMemoryEntityPydantic):
    match_id: int
    substrate_configuration: SlabStrainedSupercellConfiguration
    film_configuration: SlabStrainedSupercellConfiguration
    total_strain_percentage: Optional[float] = None
