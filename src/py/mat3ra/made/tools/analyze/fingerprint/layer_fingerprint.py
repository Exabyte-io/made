from typing import List

from pydantic import BaseModel, Field


class LayerFingerprint(BaseModel):
    min_coord: float = Field(..., description="Minimum coordinate value for the layer")
    max_coord: float = Field(..., description="Maximum coordinate value for the layer")
    elements: List[str] = Field(default_factory=list, description="Sorted unique chemical elements in the layer")

