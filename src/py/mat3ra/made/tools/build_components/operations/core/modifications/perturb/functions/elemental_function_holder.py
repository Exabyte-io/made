from typing import Any, List, Optional, Union

import sympy as sp

from mat3ra.made.periodic_table import get_atomic_mass_from_element
from .function_holder import AXIS_TO_INDEX_MAP, FunctionHolder


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

    def apply_function(
        self, coordinate: List[float], material=None, atom_index: Optional[int] = None, **kwargs: Any
    ) -> Union[float, List[float]]:
        if material is None or atom_index is None:
            raise ValueError("AtomicMassDependentFunctionHolder requires 'material' and 'atom_index' kwargs")

        element = material.basis.elements.values[atom_index]
        mass = get_atomic_mass_from_element(element)

        values = []
        for var in self.variables:
            if var == "m":
                values.append(mass)
            elif var in AXIS_TO_INDEX_MAP:
                values.append(coordinate[AXIS_TO_INDEX_MAP[var]])
            else:
                raise ValueError(f"Unknown variable: {var}")

        return self.function_numeric(*values)
