from typing import List, Callable, Optional

from mat3ra.code.entity import InMemoryEntityPydantic
import sympy as sp

from mat3ra.made.utils import AXIS_TO_INDEX_MAP


def default_function() -> sp.Expr:
    return sp.Symbol("f")


class FunctionHolder(InMemoryEntityPydantic):
    variables: List[str] = ["x"]
    symbols: List[sp.Symbol] = [sp.Symbol(var) for var in variables]
    function: sp.Expr = sp.Symbol("f")
    function_numeric: Callable = default_function
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
