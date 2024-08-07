from typing import Union, Optional

from mat3ra.made.material import Material
from .builders import (
    SlabPerturbationBuilder,
    DistancePreservingSlabPerturbationBuilder,
    CellMatchingDistancePreservingSlabPerturbationBuilder,
)
from .configuration import PerturbationConfiguration


def create_perturbation(
    configuration: PerturbationConfiguration,
    preserve_distance: Optional[bool] = False,
    builder: Union[
        SlabPerturbationBuilder,
        DistancePreservingSlabPerturbationBuilder,
        CellMatchingDistancePreservingSlabPerturbationBuilder,
        None,
    ] = None,
) -> Material:
    """
    Return a material with a perturbation applied.

    Args:
        configuration: The configuration of the perturbation to be applied.
        preserve_distance: If True, the builder that preserves the distance between atoms is used.
        builder: The builder to be used to create the perturbation.

    Returns:
        The material with the perturbation applied.
    """
    if builder is None:
        builder = SlabPerturbationBuilder()
    if preserve_distance:
        builder = CellMatchingDistancePreservingSlabPerturbationBuilder()
    return builder.get_material(configuration)
