from mat3ra.made.tools.build.defect.slab.configuration import SlabStackConfiguration


class TerraceDefectConfiguration(SlabStackConfiguration):
    cut_direction: list[int] = [1, 0, 0]
