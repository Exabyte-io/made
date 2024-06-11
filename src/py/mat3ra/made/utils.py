from typing import Any, Dict, List, Union

from mat3ra.utils.array import convert_to_array_if_not


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
