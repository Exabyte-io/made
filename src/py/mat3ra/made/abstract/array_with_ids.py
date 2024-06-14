import json


class ArrayWithIds:
    def __init__(self, array=None):
        if array is None:
            array = []
        if not isinstance(array, list):
            raise ValueError("ArrayWithIds.constructor: pass list on initialization")
        # if passed an array with ids as config, only store the values in array
        if all(isinstance(item, dict) and 'id' in item and 'value' in item for item in array):
            self.array = sorted(array, key=lambda x: x['id'])
            self.array = [element['value'] for element in self.array]
        else:
            self.array = array.copy()

    def to_json(self):
        # from ["a", "b"] to [{"id": 0, "value": "a"}, {"id": 1, "value": "b"}]
        return json.dumps([{'id': idx, 'value': el} for idx, el in enumerate(self.array)])

    def map_array_in_place(self, fn):
        if not callable(fn):
            raise ValueError("ArrayWithIds.mapArray: must pass function as argument")
        self.array = list(map(fn, self.array))

    def get_array_element_by_index(self, idx):
        return self.array[idx]

    def get_array_index_by_predicate(self, predicate):
        return next((i for i, el in enumerate(self.array) if predicate(el)), None)

    def add_element(self, el):
        value = el['value'] if isinstance(el, dict) and 'value' in el else el
        if el:
            self.array.append(value)

    def remove_element(self, el=None, idx=None):
        if idx is not None:
            _idx = idx
        else:
            _idx = self.array.index(el) if el in self.array else None
        if _idx is not None:
            self.array.pop(_idx)


