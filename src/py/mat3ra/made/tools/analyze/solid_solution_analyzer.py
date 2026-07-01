from typing import List, Optional

from mat3ra.made.tools.analyze import BaseMaterialAnalyzer
from mat3ra.made.tools.build.defective_structures.zero_dimensional.solid_solution.configuration import (
    SolidSolutionConfiguration,
)


def _most_isotropic_dimensions(total_cells: int) -> List[int]:
    best = [1, 1, total_cells]
    best_spread = total_cells - 1
    for a in range(1, total_cells + 1):
        if total_cells % a != 0:
            continue
        remaining = total_cells // a
        for b in range(a, remaining + 1):
            if remaining % b != 0:
                continue
            c = remaining // b
            if c < b:
                continue
            spread = c - a
            if spread < best_spread:
                best = [a, b, c]
                best_spread = spread
    return best


class SolidSolutionAnalyzer(BaseMaterialAnalyzer):
    """
    Plans supercell dimensions to achieve a target composition and produces a Configuration.

    Given a unit cell material and a desired substitution ratio, computes the
    smallest (most isotropic) supercell that achieves the target concentration
    within the specified tolerance.
    """

    source_element: str
    target_element: str
    target_concentration: float
    tolerance: float = 0.01
    seed: Optional[int] = None
    site_selection_method: str = "random"

    @property
    def source_element_count_per_cell(self) -> int:
        return sum(1 for el in self.material.basis.elements.values if el == self.source_element)

    @property
    def optimal_supercell_dimensions(self) -> List[int]:
        n_per_cell = self.source_element_count_per_cell
        if n_per_cell == 0:
            raise ValueError(f"No {self.source_element} atoms found in the unit cell.")

        max_cells = 512
        for total_cells in range(1, max_cells + 1):
            n_source = n_per_cell * total_cells
            n_replace = round(self.target_concentration * n_source)
            n_replace = max(0, min(n_replace, n_source))
            achieved = n_replace / n_source
            if abs(achieved - self.target_concentration) <= self.tolerance:
                return _most_isotropic_dimensions(total_cells)

        raise ValueError(
            f"Cannot achieve concentration {self.target_concentration} "
            f"within tolerance {self.tolerance} for up to {max_cells} cells."
        )

    @property
    def achievable_concentration(self) -> float:
        dims = self.optimal_supercell_dimensions
        total_cells = dims[0] * dims[1] * dims[2]
        n_source = self.source_element_count_per_cell * total_cells
        n_replace = round(self.target_concentration * n_source)
        n_replace = max(0, min(n_replace, n_source))
        return n_replace / n_source

    @property
    def number_of_substitutions(self) -> int:
        dims = self.optimal_supercell_dimensions
        total_cells = dims[0] * dims[1] * dims[2]
        n_source = self.source_element_count_per_cell * total_cells
        n_replace = round(self.target_concentration * n_source)
        return max(0, min(n_replace, n_source))

    @property
    def configuration(self) -> SolidSolutionConfiguration:
        return SolidSolutionConfiguration(
            crystal=self.material,
            source_element=self.source_element,
            target_element=self.target_element,
            concentration=self.achievable_concentration,
            supercell_dimensions=self.optimal_supercell_dimensions,
            seed=self.seed,
            site_selection_method=self.site_selection_method,
        )
