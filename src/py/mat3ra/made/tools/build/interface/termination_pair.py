from typing import Tuple
from pydantic import BaseModel


class TerminationPair(BaseModel):
    self: Tuple[str, str]

    def __init__(self, termination_pair: Tuple[str, str]):
        super().__init__(self=termination_pair)

    @property
    def film_termination(self) -> str:
        return self.self[0]

    @property
    def substrate_termination(self) -> str:
        return self.self[1]
