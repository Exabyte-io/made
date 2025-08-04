from typing import Dict, Optional

from ..enums import NanoparticleShapesEnum
from ... import BaseConfigurationPydantic


class ASEBasedNanoparticleConfiguration(BaseConfigurationPydantic):
    """
    Configuration for building a nanoparticle with bottom-up approach from atoms.

    Attributes:
        shape (NanoparticleShapes): The desired shape of the nanoparticle.
        parameters (dict): Dictionary of parameters to pass to the corresponding ASE constructor.
    """

    shape: NanoparticleShapesEnum
    parameters: Optional[Dict] = None
    element: str
    vacuum_padding: float = 10.0