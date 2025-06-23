from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.made.tools.build.slab.configuration import SlabStrainedSupercellConfiguration


class MatchedSubstrateFilmConfigurationHolder(InMemoryEntityPydantic):
    """
    Holder for matched substrate and film configurations.

    Used by interface analyzers to store paired slab configurations
    that have been processed for interface creation.
    """

    match_id: int
    substrate_configuration: SlabStrainedSupercellConfiguration
    film_configuration: SlabStrainedSupercellConfiguration
