from typing import Optional, Union

from mat3ra.made.material import Material

from .....build_components import MaterialWithBuildMetadata
from .....analyze.solid_solution_analyzer import SolidSolutionAnalyzer
from .builder import SolidSolutionBuilder


def create_solid_solution(
    material: Union[Material, MaterialWithBuildMetadata],
    source_element: str,
    target_element: str,
    concentration: float,
    seed: Optional[int] = None,
    tolerance: float = 0.01,
    site_selection_method: str = "uniform",
) -> MaterialWithBuildMetadata:
    """
    Create a solid solution by partially substituting one element for another.

    Automatically determines the optimal supercell size to achieve the target
    concentration, builds the supercell, and performs substitution.

    Args:
        material (Union[Material, MaterialWithBuildMetadata]): Unit cell crystal.
        source_element (str): Element to partially replace (e.g. "Hf").
        target_element (str): Replacement element (e.g. "Zr").
        concentration (float): Fraction of source_element to replace (0.0-1.0).
        seed (Optional[int]): Random seed for reproducible site selection.
        tolerance (float): Acceptable deviation from target concentration.
        site_selection_method (str): "random" or "uniform" (Farthest Point Sampling).

    Returns:
        MaterialWithBuildMetadata: Solid solution with full build metadata.
    """
    analyzer = SolidSolutionAnalyzer(
        material=material,
        source_element=source_element,
        target_element=target_element,
        target_concentration=concentration,
        tolerance=tolerance,
        seed=seed,
        site_selection_method=site_selection_method,
    )
    config = analyzer.configuration
    builder = SolidSolutionBuilder()
    return builder.get_material(config)
