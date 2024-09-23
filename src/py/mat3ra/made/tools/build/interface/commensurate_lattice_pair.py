import numpy as np
from pydantic import BaseModel

from .configuration import TwistedInterfaceConfiguration


class CommensurateLatticePair(BaseModel):
    class Config:
        arbitrary_types_allowed = True

    configuration: TwistedInterfaceConfiguration
    matrix1: np.ndarray
    matrix2: np.ndarray
    angle: float
    size_metric: float
