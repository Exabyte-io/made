from typing import List, Tuple, Dict, Any
import json
from .basis import Basis
from ..abstract import ArrayWithIds


class ConstrainedBasis(Basis):
    def __init__(self, config):
        super().__init__(config)
        self._constraints = ArrayWithIds(config['constraints'])

    @property
    def constraints(self):
        return self._constraints.array

    @constraints.setter
    def constraints(self, newConstraints):
        self._constraints = ArrayWithIds(newConstraints)

    def getConstraintAsArray(self):
        return self._constraints

    def getConstraintByIndex(self, idx):
        return self._constraints.getArrayElementByIndex(idx) or []

    def toJSON(self):
        return {
            **super().toJSON(),
            'constraints': self._constraints.toJSON(),
        }

    @property
    def elementsCoordinatesConstraintsArray(self):
        # Assuming _elements and other required properties/methods are defined elsewhere
        return [(element, self.getCoordinateByIndex(idx), self.getConstraintByIndex(idx), self.atomicLabelsArray[idx])
                for idx, element in enumerate(self._elements.array)]

    @property
    def atomicPositionsWithConstraints(self):
        # Assuming s.sprintf is similar to Python's format method
        return [
            "{:<4}{}".format(element + atomicLabel,
                             " ".join(["{:.9f}".format(x).strip() for x in coordinate]) + " " +
                             " ".join([str(int(x)) for x in constraint]))
            for element, coordinate, constraint, atomicLabel in self.elementsCoordinatesConstraintsArray
        ]


