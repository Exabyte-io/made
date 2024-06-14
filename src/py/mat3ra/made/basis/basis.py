from exabyte_io.periodic_table import getElectronegativity, getElementAtomicRadius
import underscore as _
import underscore.string as s
from abstract.array_with_ids import ArrayWithIds
from abstract.scalar_with_id import ObjectWithIdAndValue, ValueOrObjectArray
from constants import ATOMIC_COORD_UNITS, HASH_TOLERANCE
from lattice.lattice import Lattice, nonPeriodicLatticeScalingFactor
from lattice.types import Vector
import math
from .types import Coordinate

class Basis:
    def __init__(self, elements=["Si"], coordinates=[[0, 0, 0]], units, cell=Lattice.defaultCell, isEmpty=False, labels=[]):
        _elements = [] if isEmpty else elements
        _coordinates = [] if isEmpty else coordinates
        _units = units or Basis.unitsOptionsDefaultValue
        if not units:
            print("Basis.constructor: units are not provided => set to crystal")
        self._elements = ArrayWithIds(_elements)
        self._coordinates = ArrayWithIds(_coordinates)
        self.units = _units
        self.cell = cell
        if not _.isEmpty(labels):
            self.labels = labels

    @staticmethod
    def get_unitsOptionsConfig():
        return ATOMIC_COORD_UNITS

    @staticmethod
    def get_unitsOptionsDefaultValue():
        return ATOMIC_COORD_UNITS.crystal

    @staticmethod
    def get_defaultCell():
        return Lattice().vectorArrays

    def toJSON(self, skipRounding=False):
        json = {
            "elements": self.elements,
            "coordinates": self.coordinates if skipRounding else self.coordinatesRounded,
            "units": self.units,
            "cell": self.cell if skipRounding else self.cellRounded
        }
        if not _.isEmpty(self.labels):
            return JSON.parse(
                JSON.stringify({
                    **json,
                    "labels": self.labels
                })
            )
        return JSON.parse(JSON.stringify(json))

    @property
    def coordinatesRounded(self):
        return [
            {
                "id": coordinate["id"],
                "value": [math.precise(math.roundToZero(x)) for x in coordinate["value"]]
            }
            for coordinate in self.coordinates
        ]

    @property
    def cellRounded(self):
        return [
            [math.precise(math.roundToZero(x)) for x in vector]
            for vector in self.cell
        ]

    def clone(self, extraContext=None):
        return Basis(
            elements=self.toJSON()["elements"],
            coordinates=self.toJSON()["coordinates"],
            units=self.units,
            cell=self.cell,
            isEmpty=False,
            labels=self.labels
        )

    def getElementByIndex(self, idx):
        return self._elements.getArrayElementByIndex(idx)

    def getCoordinateByIndex(self, idx):
        return self._coordinates.getArrayElementByIndex(idx)

    @property
    def elementsArray(self):
        return self._elements.array

    @property
    def elements(self):
        return self._elements.toJSON()

    @elements.setter
    def elements(self, elementsArray):
        self._elements = ArrayWithIds(elementsArray)

    def getElementsAsObject(self):
        return self._elements.toJSON()

    @property
    def coordinates(self):
        return self._coordinates.toJSON()

    @coordinates.setter
    def coordinates(self, coordinatesNestedArray):
        self._coordinates = ArrayWithIds(coordinatesNestedArray)

    @property
    def coordinatesAsArray(self):
        return self._coordinates.array

    @property
    def isInCrystalUnits(self):
        return self.units == ATOMIC_COORD_UNITS.crystal

    @property
    def isInCartesianUnits(self):
        return self.units == ATOMIC_COORD_UNITS.cartesian

    def toCartesian(self):
        unitCell = self.cell
        if self.units == ATOMIC_COORD_UNITS.cartesian:
            return
        self._coordinates.mapArrayInPlace(
            lambda point: math.multiply(point, unitCell) as Coordinate
        )
        self.units = ATOMIC_COORD_UNITS.cartesian

    def toCrystal(self):
        unitCell = self.cell
        if self.units == ATOMIC_COORD_UNITS.crystal:
            return
        self._coordinates.mapArrayInPlace(
            lambda point: math.multiply(point, math.inv(unitCell)) as Coordinate
        )
        self.units = ATOMIC_COORD_UNITS.crystal

    def toStandardRepresentation(self):
        self.toCrystal()
        self._coordinates.mapArrayInPlace(
            lambda point: [math.mod(x) for x in point]
        )

    @property
    def standardRepresentation(self):
        originalUnits = self.units
        self.toStandardRepresentation()
        result = self.toJSON()
        if originalUnits != ATOMIC_COORD_UNITS.crystal:
            self.toCartesian()
        return result

    def addAtom(self, element="Si", coordinate=[0.5, 0.5, 0.5]):
        self._elements.addElement(element)
        self._coordinates.addElement(coordinate)

    def removeAtom(self, element=None, coordinate=None, id=None):
        if element and len(coordinate) > 0:
            self._elements.removeElement(element, id)
            self._coordinates.removeElement(coordinate, id)

    @property
    def uniqueElements(self):
        return _.unique(self._elements.array)

    @property
    def uniqueElementCountsSortedByElectronegativity(self):
        return _.chain(self.elements) \
            .sortBy("value") \
            .sortBy(lambda x: getElectronegativity(x["value"])) \
            .countBy(lambda element: element["value"]) \
            .value()

    @property
    def elementCounts(self):
        elementCounts = []
        elements = self.getElementsAsObject()
        for index, element in enumerate(elements):
            previousElement = elements[index - 1]
            if previousElement and previousElement["value"] == element["value"]:
                previousElementCount = elementCounts[-1]
                previousElementCount["count"] += 1
            else:
                elementCounts.append({
                    "count": 1,
                    "value": element["value"]
                })
        return elementCounts

    @property
    def formula(self):
        counts = self.uniqueElementCountsSortedByElectronegativity
        countsValues = _.values(counts)
        gcd = math.gcd(*countsValues) if len(countsValues) > 1 else countsValues[0]
        return "".join([
            f"{x[0]}{'' if x[1] / gcd == 1 else x[1] / gcd}"
            for x in _.pairs(counts)
        ])

    @property
    def unitCellFormula(self):
        counts = self.uniqueElementCountsSortedByElectronegativity
        return "".join([
            f"{x[0]}{'' if x[1] == 1 else x[1]}"
            for x in _.pairs(counts)
        ])

    @property
    def elementsAndCoordinatesArray(self):
        return [
            [element, self.getCoordinateByIndex(idx)]
            for idx, element in enumerate(self.elementsArray)
        ]

    @property
    def elementsAndCoordinatesAndLabelsArray(self):
        return [
            [element, self.getCoordinateByIndex(idx), self.atomicLabelsArray[idx]]
            for idx, element in enumerate(self.elementsArray)
        ]

    @property
    def atomicLabelsArray(self):
        labelsArray = [""] * len(self.elements)
        if self.labels:
            for item in self.labels:
                labelsArray[item["id"]] = str(item["value"])
        return labelsArray

    @property
    def elementsWithLabelsArray(self):
        return [
            f"{symbol}{label}"
            for symbol, label in zip(self.elementsArray, self.atomicLabelsArray)
        ]

    @property
    def atomicPositions(self):
        return [
            f"{element}{label} {' '.join([s.sprintf('%14.9f', x).strip() for x in coordinate])}"
            for element, coordinate, label in self.elementsAndCoordinatesAndLabelsArray
        ]

    @property
    def nAtoms(self):
        return len(self._elements.array)

    def isEqualTo(self, anotherBasisClsInstance):
        return self.hashString == anotherBasisClsInstance.hashString

    def hasEquivalentCellTo(self, anotherBasisClsInstance):
        return all([
            math.vEqualWithTolerance(vector, anotherBasisClsInstance.cell[idx])
            for idx, vector in enumerate(self.cell)
        ])

    def getMinimumLatticeSize(self, latticeScalingFactor=nonPeriodicLatticeScalingFactor):
        latticeSizeAdditiveContribution = 0
        if len(self._elements.array) == 1:
            elementSymbol = self._elements.getArrayElementByIndex(0)
            latticeSizeAdditiveContribution = getElementAtomicRadius(elementSymbol)
        moleculeLatticeSize = self.maxPairwiseDistance * latticeScalingFactor
        latticeSize = latticeSizeAdditiveContribution + moleculeLatticeSize
        return round(latticeSize, 4)

    def getOverlappingAtoms(self):
        self.toCartesian()
        coordinates = self.coordinates
        elements = self.elements
        overlaps = []
        overlapCoefficient = 0.75
        for i, entry1 in enumerate(coordinates):
            for j in range(i + 1, len(coordinates)):
                entry2 = coordinates[j]
                el1 = elements[i]["value"]
                el2 = elements[j]["value"]
                tolerance = overlapCoefficient * (getElementAtomicRadius(el1) + getElementAtomicRadius(el2))
                distance = math.vDist(entry1["value"], entry2["value"])
                if distance < tolerance:
                    overlaps.append({
                        "id1": i,
                        "id2": j,
                        "element1": el1,
                        "element2": el2
                    })
        self.toCrystal()
        return overlaps

    @property
    def maxPairwiseDistance(self):
        originalUnits = self.units
        self.toCartesian()
        maxDistance = 0
        if len(self._elements.array) >= 2:
            for i in range(len(self._elements.array)):
                for j in range(i + 1, len(self._elements.array)):
                    distance = math.vDist(
                        self._coordinates.getArrayElementByIndex(i),
                        self._coordinates.getArrayElementByIndex(j)
                    )
                    if distance > maxDistance:
                        maxDistance = distance
        if originalUnits:
            self.toCrystal()
        return maxDistance

    @property
    def hashString(self):
        originallyInCartesianUnits = self.isInCartesianUnits
        self.toCrystal()
        hashString = self.getAsSortedString()
        if originallyInCartesianUnits:
            self.toCartesian()
        return hashString

    def getAsSortedString(self):
        clsInstance = Basis(self.toJSON())
        clsInstance.toStandardRepresentation()
        standardRep = [
            f"{element}{atomicLabel} {' '.join([str(math.round(x, HASH_TOLERANCE)) for x in coordinate])}"
            for element, coordinate, atomicLabel in clsInstance.elementsAndCoordinatesAndLabelsArray
        ]
        return ";".join(sorted(standardRep)) + ";"


