from typing import List, Optional

from ...build_components.entities.auxiliary.two_dimensional.termination import Termination


def select_slab_termination(terminations: List[Termination], formula: Optional[str] = None) -> Termination:
    if not terminations:
        raise ValueError("No terminations available.")
    if formula is None:
        return terminations[0]
    for termination in terminations:
        if termination.formula == formula:
            return termination
    raise ValueError(f"Termination with formula {formula} not found in available terminations: {terminations}")
