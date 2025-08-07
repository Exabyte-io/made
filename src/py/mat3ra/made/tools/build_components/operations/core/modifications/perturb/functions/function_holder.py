from typing import Any, Callable, List, Optional, Union

import sympy as sp
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.made.utils import AXIS_TO_INDEX_MAP
from pydantic import field_serializer
from sympy import parse_expr


def default_function() -> sp.Expr:
    return sp.Symbol("f")


class FunctionHolder(InMemoryEntityPydantic):
    variables: List[str] = ["x"]
    symbols: List[sp.Symbol] = sp.symbols(variables)
    function: sp.Expr = sp.Symbol("f")
    function_numeric: Callable = sp.lambdify(sp.symbols(variables), sp.Symbol("f"), modules=["numpy"])
    derivatives_numeric: dict = {}

    def __init__(self, function: Union[sp.Expr, str], variables: Optional[List[str]] = None, **data: Any):
        # normalize string → Expr
        expr = self._to_expr(function)

        # if no variables list given, infer from free symbols (defaults to x,y,z)
        if variables is None:
            vs = sorted(expr.free_symbols, key=lambda s: s.name)
            variables = [str(v) for v in vs] or ["x", "y", "z"]
        super().__init__(**data)

        self.variables = variables
        self.symbols = sp.symbols(variables)
        self.function = expr

        self.function_numeric = sp.lambdify(self.symbols, self.function, modules=["numpy"])
        self.derivatives_numeric = {
            var: sp.lambdify(sp.symbols(variables), sp.diff(self.function, var), modules=["numpy"])
            for var in self.variables
        }

    @staticmethod
    def _to_expr(expr_or_str: Union[sp.Expr, str]) -> sp.Expr:
        if isinstance(expr_or_str, sp.Expr):
            return expr_or_str
        if isinstance(expr_or_str, str):
            # We can swap in sympy.sympify for slightly laxer parsing,
            # or configure parse_expr with custom transformations.
            return parse_expr(expr_or_str, evaluate=True)
        raise TypeError(f"Expected sympy.Expr or str, got {type(expr_or_str)}")

    @property
    def function_str(self) -> str:
        return str(self.function)

    def apply_function(self, coordinate: List[float]) -> float:
        values = [coordinate[AXIS_TO_INDEX_MAP[var]] for var in self.variables]
        return self.function_numeric(*values)

    def calculate_derivative(self, coordinate: List[float], axis: str) -> float:
        if axis in self.variables:
            values = [coordinate[AXIS_TO_INDEX_MAP[var]] for var in self.variables]
            return self.derivatives_numeric[axis](*values)
        else:
            return 0

    @field_serializer("function")
    def serialize_function(self, value: sp.Expr) -> str:
        return str(value)

    @field_serializer("symbols")
    def serialize_symbols(self, value: List[sp.Symbol]) -> List[str]:
        return [str(symbol) for symbol in value]

    @field_serializer("function_numeric")
    def serialize_function_numeric(self, value: Callable) -> str:
        return "lambdified_function"

    @field_serializer("derivatives_numeric")
    def serialize_derivatives_numeric(self, value: dict) -> dict:
        return {k: "lambdified_derivative" for k, v in value.items()}
