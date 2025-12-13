from typing import Any, Callable, List, Optional, Union

import sympy as sp
from mat3ra.code.entity import InMemoryEntityPydantic
from pydantic import field_serializer
from sympy import parse_expr

from mat3ra.made.utils import AXIS_TO_INDEX_MAP


class ElementalFunctionHolder(InMemoryEntityPydantic):
    variables: List[str] = ["x", "y", "z", "m"]
    symbols: List[sp.Symbol] = sp.symbols(["x", "y", "z", "m"])
    function: sp.Expr = sp.Symbol("f")
    function_numeric: Callable = None
    atom_masses: List[float] = []

    def __init__(
        self,
        function: Union[sp.Expr, str],
        atom_masses: List[float],
        variables: Optional[List[str]] = None,
        **data: Any,
    ):
        expr = self._to_expr(function)

        if variables is None:
            vs = sorted(expr.free_symbols, key=lambda s: s.name)
            variables = [str(v) for v in vs] or ["x", "y", "z", "m"]

        super().__init__(**data)

        self.variables = variables
        self.symbols = sp.symbols(variables)
        self.function = expr
        self.atom_masses = atom_masses

        self.function_numeric = sp.lambdify(self.symbols, self.function, modules=["numpy"])

    @staticmethod
    def _to_expr(expr_or_str: Union[sp.Expr, str]) -> sp.Expr:
        if isinstance(expr_or_str, sp.Expr):
            return expr_or_str
        if isinstance(expr_or_str, str):
            return parse_expr(expr_or_str, evaluate=True)
        raise TypeError(f"Expected sympy.Expr or str, got {type(expr_or_str)}")

    @property
    def function_str(self) -> str:
        return str(self.function)

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

    @field_serializer("function")
    def serialize_function(self, value: sp.Expr) -> str:
        return str(value)

    @field_serializer("symbols")
    def serialize_symbols(self, value: List[sp.Symbol]) -> List[str]:
        return [str(symbol) for symbol in value]

    @field_serializer("function_numeric")
    def serialize_function_numeric(self, value: Callable) -> str:
        return "lambdified_function"
