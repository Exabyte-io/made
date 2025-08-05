from typing import Dict, List

from pydantic import BaseModel


class CoordinateCondition(BaseModel):
    class Config:
        arbitrary_types_allowed = True

    def condition(self, coordinate: List[float]) -> bool:
        raise NotImplementedError

    def get_json(self) -> Dict:
        json = {"type": self.__class__.__name__}
        json.update(self.model_dump())
        return json
