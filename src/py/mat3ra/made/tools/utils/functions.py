from typing import List, Callable, Optional, Any

import sympy as sp
from mat3ra.code.entity import InMemoryEntityPydantic
from pydantic import field_serializer

from mat3ra.made.utils import AXIS_TO_INDEX_MAP


def default_function() -> sp.Expr:
    return sp.Symbol("f")


class FunctionHolder(InMemoryEntityPydantic):
    variables: List[str] = ["x"]
    symbols: List[sp.Symbol] = sp.symbols(variables)
    function: sp.Expr = sp.Symbol("f")
    function_numeric: Callable = sp.lambdify(sp.symbols(variables), sp.Symbol("f"), modules=["numpy"])
    derivatives_numeric: dict = {}

    def __init__(self, function: Optional[sp.Expr] = None, variables: Optional[List[str]] = None, **data: Any):
        """
        Initializes with a function involving multiple variables.
        """
        if function is None:
            function = default_function()
        if variables is None:
            variables = ["x", "y", "z"]
        super().__init__(**data)
        self.variables = variables
        self.symbols = sp.symbols(variables)
        self.function = function
        self.function_numeric = sp.lambdify(self.symbols, self.function, modules=["numpy"])
        self.derivatives_numeric = {
            var: sp.lambdify(self.symbols, sp.diff(self.function, var), modules=["numpy"]) for var in variables
        }

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
