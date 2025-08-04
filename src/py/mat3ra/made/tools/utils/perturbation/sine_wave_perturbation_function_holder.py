from typing import Any

import sympy as sp

from .perturbation_function_holder import PerturbationFunctionHolder


class SineWavePerturbationFunctionHolder(PerturbationFunctionHolder):
    amplitude: float = 0.05
    wavelength: float = 1
    phase: float = 0
    axis: str = "x"

    def __init__(
        self,
        amplitude: float = 0.05,
        wavelength: float = 1,
        phase: float = 0,
        axis: str = "x",
        **data: Any,
    ):
        w = sp.Symbol(axis)
        function = amplitude * sp.sin(2 * sp.pi * w / wavelength + phase)
        variables = [axis]
        super().__init__(function=function, variables=variables, **data)
        self.amplitude = amplitude
        self.wavelength = wavelength
        self.phase = phase
        self.axis = axis
