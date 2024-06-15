import json
from typing import Any, Dict, List, Union, Callable

from mat3ra.utils.array import convert_to_array_if_not
from mat3ra.utils.mixins import RoundNumericValuesMixin
from pydantic import BaseModel


# TODO: move to a more general location
def map_array_to_array_with_id_value(array: List[Any], remove_none: bool = False) -> List[Any]:
    full_array = [{"id": i, "value": item} for i, item in enumerate(array)]
    if remove_none:
        return list(filter(lambda x: x["value"] is not None, full_array))
    return full_array


def map_array_with_id_value_to_array(array: List[Dict[str, Any]]) -> List[Any]:
    return [item["value"] for item in array]


def get_array_with_id_value_element_value_by_index(array: List[Dict[str, Any]], index: int = 0) -> List[Any]:
    return map_array_with_id_value_to_array(array)[index]


def filter_array_with_id_value_by_values(
    array: List[Dict[str, Any]], values: Union[List[Any], Any]
) -> List[Dict[str, Any]]:
    values = convert_to_array_if_not(values)
    return [item for item in array if item["value"] in values]
    # Alternative implementation:
    # return list(filter(lambda x: x["value"] in values, array))


def filter_array_with_id_value_by_ids(
    array: List[Dict[str, Any]], ids: Union[List[int], List[str], int, str]
) -> List[Dict[str, Any]]:
    int_ids = list(map(lambda i: int(i), convert_to_array_if_not(ids)))
    return [item for item in array if item["id"] in int_ids]
    # Alternative implementation:
    # return list(filter(lambda x: x["id"] in ids, array))


def are_arrays_equal_by_id_value(array1: List[Dict[str, Any]], array2: List[Dict[str, Any]]) -> bool:
    return map_array_with_id_value_to_array(array1) == map_array_with_id_value_to_array(array2)


class ValueWithId(BaseModel):
    id: int = 0
    value: Any = None

    def to_dict(self) -> Dict[str, Any]:
        return {"id": self.id, "value": self.value}

    def to_json(self) -> str:
        return json.dumps(self.to_dict())

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ValueWithId):
            return False
        return self.value == other.value and self.id == other.id


class ArrayWithIds(RoundNumericValuesMixin, BaseModel):
    array: List[Any] = []

    def to_array_of_values_with_ids(self) -> List[ValueWithId]:
        return [ValueWithId(id=index, value=item) for index, item in enumerate(self.array)]

    def get_element_value_by_index(self, index: int = 0) -> List[Any]:
        return self.array[index]

    def filter_by_values(self, values: Union[List[Any], Any]):
        values = convert_to_array_if_not(values)
        self.array = [item for item in self.array if item in values]

    def filter_by_ids(self, ids: Union[List[int], List[str], int, str]):
        int_ids = list(map(lambda i: int(i), convert_to_array_if_not(ids)))
        self.array = [item for index, item in enumerate(self.array) if index in int_ids]

    def to_json(self, skip_rounding=True) -> str:
        array_with_ids_values = [
            self.round_array_or_number(item) if skip_rounding else item for item in self.to_array_of_values_with_ids()
        ]
        return json.dumps([item.to_json() for item in array_with_ids_values])

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ArrayWithIds):
            return False
        return self.array == other.array

    def map_array_in_place(self, func: Callable):
        self.array = list(map(func, self.array))

    def add_item(self, element: Any):
        self.array.append(element)

    def remove_item(self, index: int = 0):
        self.array.pop(index)
