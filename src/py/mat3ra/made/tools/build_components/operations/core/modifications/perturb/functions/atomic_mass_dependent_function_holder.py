from typing import Any, List, Optional, Union

import sympy as sp
from mat3ra.periodic_table.helpers import get_atomic_mass_from_element

from .function_holder import FunctionHolder


class AtomicMassDependentFunctionHolder(FunctionHolder):
    variables: List[str] = ["x", "y", "z", "m"]

    def __init__(
        self,
        function: Union[sp.Expr, str],
        variables: Optional[List[str]] = None,
        **data: Any,
    ):
        if variables is None:
            expr = self._to_expr(function)
            vs = sorted(expr.free_symbols, key=lambda s: s.name)
            variables = [str(v) for v in vs] or ["x", "y", "z", "m"]

        super().__init__(function=function, variables=variables, **data)

    @staticmethod
    def get_atomic_mass(coordinate: List[float], material) -> float:
        if material is None:
            raise ValueError("Material is required to extract atomic mass")

        atom_id = material.basis.coordinates.get_element_id_by_value(coordinate)
        element = material.basis.elements.get_element_value_by_index(atom_id)
        return get_atomic_mass_from_element(element)
