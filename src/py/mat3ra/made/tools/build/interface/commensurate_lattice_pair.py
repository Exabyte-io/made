import numpy as np
from pydantic import BaseModel

from .configuration import TwistedInterfaceConfiguration


class CommensurateLatticePair(BaseModel):
    """
    Commensurate lattice pair model.

    Attributes:
        configuration (TwistedInterfaceConfiguration): The configuration of the twisted interface.
        matrix1 (np.ndarray): The supercell 2D matrix for the first lattice.
        matrix2 (np.ndarray): The supercell 2D matrix for the second lattice.
        angle (float): The angle between the two lattices, in degrees.
        size_metric (float): The size metric of the resulting supercell, in arbitrary units.
    """

    class Config:
        arbitrary_types_allowed = True

    configuration: TwistedInterfaceConfiguration
    matrix1: np.ndarray
    matrix2: np.ndarray
    angle: float
    size_metric: float
