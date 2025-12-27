from typing import Any, List, Optional

import numpy as np
import sympy as sp
from pydantic import Field

from .atomic_mass_dependent_function_holder import AtomicMassDependentFunctionHolder

DEFAULT_CONVERSION_CONSTANT = 2e-3


class MaxwellBoltzmannDisplacementHolder(AtomicMassDependentFunctionHolder):
    disorder_parameter: float = Field(exclude=True)
    random_seed: Optional[int] = Field(default=None, exclude=True)
    random_state: Any = Field(default=None, exclude=True)
    is_mass_used: bool = Field(default=True, exclude=True)
    conversion_constant: float = Field(default=DEFAULT_CONVERSION_CONSTANT, exclude=True)

    def __init__(
        self,
        disorder_parameter: float,
        random_seed: Optional[int] = None,
        is_mass_used: bool = True,
        conversion_constant: float = DEFAULT_CONVERSION_CONSTANT,
    ):
        calibrated_disorder_parameter = disorder_parameter * conversion_constant
        random_state = np.random.RandomState(random_seed) if random_seed is not None else np.random
        function_expr = sp.Symbol("f")

        super().__init__(
            function=function_expr,
            disorder_parameter=calibrated_disorder_parameter,
            random_state=random_state,
            is_mass_used=is_mass_used,
            conversion_constant=conversion_constant,
        )

    def apply_function(self, coordinate, material=None) -> List[float]:
        if material is None:
            raise ValueError("MaxwellBoltzmannDisplacementHolder requires 'material' kwargs")

        if self.is_mass_used:
            mass = self.get_atomic_mass(coordinate, material)
            variance = self.disorder_parameter / mass
        else:
            variance = self.disorder_parameter

        std_dev = np.sqrt(variance)
        displacement = self.random_state.normal(0.0, std_dev, size=3)
        return displacement.tolist()
