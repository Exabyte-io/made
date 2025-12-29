from typing import Any, List, Optional

import numpy as np
from pydantic import Field, model_validator

from .atomic_mass_dependent_function_holder import AtomicMassDependentFunctionHolder

DEFAULT_DISORDER_PARAMETER = 1.0


class MaxwellBoltzmannDisplacementHolder(AtomicMassDependentFunctionHolder):
    disorder_parameter: float = Field(
        default=DEFAULT_DISORDER_PARAMETER,
        exclude=True,
        description="Disorder parameter. Can be viewed as effective temperature in eV.",
    )
    random_seed: Optional[int] = Field(default=None, exclude=True)
    random_state: Any = Field(default=None, exclude=True)
    is_mass_used: bool = Field(default=True, exclude=True)

    @model_validator(mode="after")
    def setup_random_state(self):
        if self.random_state is None:
            self.random_state = np.random.RandomState(self.random_seed) if self.random_seed is not None else np.random
        return self

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
