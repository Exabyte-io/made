from typing import Optional

from mat3ra.code.entity import InMemoryEntityPydantic

from ....build.pristine_structures.two_dimensional.slab_strained_supercell.configuration import (
    SlabStrainedSupercellConfiguration,
)


class MatchedSubstrateFilmConfigurationHolder(InMemoryEntityPydantic):
    match_id: int
    substrate_configuration: SlabStrainedSupercellConfiguration
    film_configuration: SlabStrainedSupercellConfiguration
    total_strain_percentage: Optional[float] = None
