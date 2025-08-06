from typing import Optional

from mat3ra.code.entity import InMemoryEntityPydantic

from ..build_components.entities.auxiliary.two_dimensional.termination import Termination


class TerminationHolder(InMemoryEntityPydantic):
    termination_with_vacuum: Termination
    termination_without_vacuum: Optional[Termination]
    shift_with_vacuum: float
    shift_without_vacuum: Optional[float]
