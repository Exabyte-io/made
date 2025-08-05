from typing import Callable

from pydantic import BaseModel

from ..interaction_functions import sum_of_inverse_distances_squared


class MaterialCalculatorParameters(BaseModel):
    """
    Defines the parameters for a material calculator.

    Args:
        interaction_function (Callable): A function used to calculate the interaction metric between
        sets of coordinates. The default function is sum_of_inverse_distances_squared.
    """

    interaction_function: Callable = sum_of_inverse_distances_squared
