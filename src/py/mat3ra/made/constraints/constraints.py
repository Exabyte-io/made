from ..abstract.array_with_ids import ArrayWithIds


class AtomicConstraints:
    def __init__(self, values=None):
        self.name = "atomic_constraints"
        self.values = ArrayWithIds(values or [])

    @staticmethod
    def fromArray(array):
        return AtomicConstraints(values=array)

    def toJSON(self):
        return {
            "name": self.name,
            "values": self.values.toJSON(),
        }

    def getByIndex(self, idx):
        return self.values.getArrayElementByIndex(idx) or []

    def getAsStringByIndex(self, idx, mapFn=None):
        if mapFn is None:
            mapFn = lambda val: "1" if val else "0"
        return " ".join(map(mapFn, self.getByIndex(idx)))


