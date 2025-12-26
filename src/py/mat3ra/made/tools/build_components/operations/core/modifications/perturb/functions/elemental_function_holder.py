from typing import Any, List, Optional, Union

import sympy as sp

from .function_holder import AXIS_TO_INDEX_MAP, FunctionHolder


class AtomicMassDependentFunctionHolder(FunctionHolder):
    variables: List[str] = ["x", "y", "z", "m"]
    atom_masses: List[float] = []

    def __init__(
        self,
        function: Union[sp.Expr, str],
        atom_masses: List[float],
        variables: Optional[List[str]] = None,
        **data: Any,
    ):
        if variables is None:
            expr = self._to_expr(function)
            vs = sorted(expr.free_symbols, key=lambda s: s.name)
            variables = [str(v) for v in vs] or ["x", "y", "z", "m"]

        super().__init__(function=function, variables=variables, atom_masses=atom_masses, **data)

    def apply_function(self, coordinate: List[float], atom_index: int) -> Union[float, List[float]]:
        if atom_index < 0 or atom_index >= len(self.atom_masses):
            raise ValueError(f"Atom index {atom_index} out of range [0, {len(self.atom_masses)})")

        mass = self.atom_masses[atom_index]
        values = []
        for var in self.variables:
            if var == "m":
                values.append(mass)
            elif var in AXIS_TO_INDEX_MAP:
                values.append(coordinate[AXIS_TO_INDEX_MAP[var]])
            else:
                raise ValueError(f"Unknown variable: {var}")

        return self.function_numeric(*values)
