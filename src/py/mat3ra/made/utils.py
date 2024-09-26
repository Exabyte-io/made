import json
from typing import Any, Callable, Dict, List, Optional, Union

import numpy as np
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


def get_center_of_coordinates(coordinates: List[List[float]]) -> List[float]:
    """
    Calculate the center of the coordinates.

    Args:
        coordinates (List[List[float]]): The list of coordinates.

    Returns:
        List[float]: The center of the coordinates.
    """
    return np.mean(np.array(coordinates), axis=0).tolist()


def get_overlapping_coordinates(
    coordinate: List[float],
    coordinates: List[List[float]],
    threshold: float = 0.01,
) -> List[List[float]]:
    """
    Find coordinates that are within a certain threshold of a given coordinate.

    Args:
        coordinate (List[float]): The coordinate.
        coordinates (List[List[float]]): The list of coordinates.
        threshold (float): The threshold for the distance, in the units of the coordinates.

    Returns:
        List[List[float]]: The list of overlapping coordinates.
    """
    return [c for c in coordinates if np.linalg.norm(np.array(c) - np.array(coordinate)) < threshold]


def create_2d_supercell_matrices(max_search: int) -> List[np.ndarray]:
    """
    Create a list of 2D supercell matrices within a maximum search range.

    Filtering conditions:
    - Non-zero area constraint
    - Positive determinant (to exclude mirroring transformations)

    Args:
        max_search: The maximum search range.
    Returns:
        List[np.ndarray]: The list of supercell matrices.
    """
    matrices = []
    for s11 in range(-max_search, max_search + 1):
        for s12 in range(-max_search, max_search + 1):
            for s21 in range(-max_search, max_search + 1):
                for s22 in range(-max_search, max_search + 1):
                    matrix = np.array([[s11, s12], [s21, s22]])
                    determinant = np.linalg.det(matrix)
                    if determinant == 0 or determinant < 0:
                        continue
                    matrices.append(matrix)
    return matrices


def get_angle_from_rotation_matrix_2d(
    matrix: np.ndarray, zero_tolerance: float = 1e-6, round_digits: int = 3
) -> Union[float, None]:
    """
    Get the angle from a 2x2 rotation matrix in degrees if it's a pure rotation matrix.
    Args:
        matrix: The 2x2 rotation matrix.
        zero_tolerance: The zero tolerance for the determinant.
        round_digits: The number of digits to round the angle.

    Returns:
        Union[float, None]: The angle in degrees if it's a pure rotation matrix, otherwise None.
    """
    if matrix.shape != (2, 2):
        return None
    if np.abs(np.linalg.det(matrix) - 1) > zero_tolerance:
        return None
    if not np.all(np.abs(matrix) <= 1):
        return None
    # Check if it's in form of rotation matrix [cos(theta), -sin(theta); sin(theta), cos(theta)]
    if not np.allclose(matrix @ matrix.T, np.eye(2), atol=zero_tolerance):
        return None
    cos_theta = matrix[0, 0]
    sin_theta = matrix[1, 0]
    angle_rad = np.arctan2(sin_theta, cos_theta)
    angle_deg = np.round(np.degrees(angle_rad), round_digits)
    return angle_deg


class ValueWithId(RoundNumericValuesMixin, BaseModel):
    id: int = 0
    value: Any = None

    def to_dict(self, skip_rounding: bool = True) -> Dict[str, Any]:
        rounded_value = self.round_array_or_number(self.value) if not skip_rounding else self.value
        return {"id": self.id, "value": rounded_value}

    def to_json(self, skip_rounding: bool = True) -> str:
        return json.dumps(self.to_dict(skip_rounding=skip_rounding))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ValueWithId):
            return False
        return self.value == other.value and self.id == other.id


class ArrayWithIds(RoundNumericValuesMixin, BaseModel):
    values: List[Any] = []
    ids: List[int] = []

    @classmethod
    def from_values(cls, values: List[Any]) -> "ArrayWithIds":
        try:
            ids = list(range(len(values)))
            return cls(values=values, ids=ids)
        except KeyError:
            raise ValueError("Values must be a list")

    @classmethod
    def from_list_of_dicts(cls, list_of_dicts: List[Dict[str, Any]]) -> "ArrayWithIds":
        try:
            values = [item["value"] for item in list_of_dicts]
            ids = [item["id"] for item in list_of_dicts]
            return cls(values=values, ids=ids)
        except KeyError:
            raise ValueError("List of dictionaries must contain 'id' and 'value' keys")

    def to_array_of_values_with_ids(self) -> List[ValueWithId]:
        return [ValueWithId(id=id, value=item) for id, item in zip(self.ids, self.values)]

    def get_element_value_by_index(self, index: int) -> Any:
        return self.values[index] if index < len(self.values) else None

    def get_element_id_by_value(self, value: Any) -> Optional[int]:
        try:
            return self.ids[self.values.index(value)]
        except ValueError:
            return None

    def filter_by_values(self, values: Union[List[Any], Any]):
        def make_hashable(value):
            return tuple(value) if isinstance(value, list) else value

        values_to_keep = set(make_hashable(v) for v in values) if isinstance(values, list) else {make_hashable(values)}
        filtered_items = [(v, i) for v, i in zip(self.values, self.ids) if make_hashable(v) in values_to_keep]
        if filtered_items:
            values_unpacked, ids_unpacked = zip(*filtered_items)
            self.values = list(values_unpacked)
            self.ids = list(ids_unpacked)
        else:
            self.values = []
            self.ids = []

    def filter_by_indices(self, indices: Union[List[int], int]):
        index_set = set(indices) if isinstance(indices, list) else {indices}
        self.values = [self.values[i] for i in range(len(self.values)) if i in index_set]
        self.ids = [self.ids[i] for i in range(len(self.ids)) if i in index_set]

    def filter_by_ids(self, ids: Union[List[int], int]):
        if isinstance(ids, int):
            ids = [ids]
        ids_set = set(ids)
        keep_indices = [index for index, id_ in enumerate(self.ids) if id_ in ids_set]
        self.values = [self.values[index] for index in keep_indices]
        self.ids = [self.ids[index] for index in keep_indices]

    def to_dict(self) -> List[Dict[str, Any]]:
        return [
            item.to_dict() if isinstance(item, ValueWithId) else ValueWithId(id=index, value=item).to_dict()
            for index, item in enumerate(self.values)
        ]

    def to_json(self, skip_rounding=True) -> str:
        array_with_ids_values = [
            ValueWithId(id=item.id, value=self.round_array_or_number(item.value)) if not skip_rounding else item
            for item in self.to_array_of_values_with_ids()
        ]
        return json.loads(json.dumps([item.to_dict() for item in array_with_ids_values]))

    def __eq__(self, other: object) -> bool:
        return isinstance(other, ArrayWithIds) and self.values == other.values and self.ids == other.ids

    def map_array_in_place(self, func: Callable):
        self.values = list(map(func, self.values))

    def add_item(self, element: Any, id: Optional[int] = None):
        if id is None:
            new_id = max(self.ids, default=-1) + 1
        else:
            new_id = id
        self.values.append(element)
        self.ids.append(new_id)

    def remove_item(self, index: int, id: Optional[int] = None):
        if id is not None:
            try:
                index = self.ids.index(id)
            except ValueError:
                raise ValueError("ID not found in the list")
        if index < len(self.values):
            del self.values[index]
            del self.ids[index]
        else:
            raise IndexError("Index out of range")
