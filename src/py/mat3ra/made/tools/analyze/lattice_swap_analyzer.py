from typing import List, Tuple, Union, Optional

import numpy as np
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.tools.analyze.rotation_analyzer import MaterialRotationAnalyzer
from pydantic import BaseModel, Field

from mat3ra.made.material import Material
from .basis.analyzer import BasisMaterialAnalyzer
from .fingerprint import MaterialFingerprintAllAxes
from ..build_components.metadata.material_with_build_metadata import MaterialWithBuildMetadata
from ...lattice import Lattice


class LatticeSwapParameters(BaseModel):
    model_config = {"arbitrary_types_allowed": True}

    permutation: List[Tuple[int, int]] = Field(..., description="List of (source_idx, sign) tuples for each new vector")


class LatticeSwapDetectionResult(LatticeSwapParameters):
    is_swapped: bool = Field(..., description="Whether a lattice swap was detected")
    permutation: List[Tuple[int, int]] = Field(..., description="Permutation that was detected")
    new_lattice: Lattice = Field(..., description="New lattice configuration with swapped parameters")
    confidence: float = Field(..., description="Confidence score (0.0 to 1.0)")


class MaterialLatticeSwapAnalyzer(MaterialRotationAnalyzer):
    """
    Analyzer to detect lattice vector swaps/permutations between materials.

    This detects when lattice vectors have been reoriented (e.g., a->a, b->c, c->-b)
    rather than the basis being rotated within the lattice.
    """

    material: Union[Material, MaterialWithBuildMetadata]
    tolerance: float = 0.01

    # Let's reuse the rotation detection logic for fingerprints, then when we get a rotation we convert it to a corresponding swap of lattice vectors and angles
