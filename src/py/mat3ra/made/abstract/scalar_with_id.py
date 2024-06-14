from typing import Union, List, TypeVar, Generic, Dict

T = TypeVar('T')


class ObjectWithIdAndValue(Generic[T]):
    def __init__(self, id: int, value: T):
        self.id = id
        self.value = value


ValueOrObject = Union[ObjectWithIdAndValue[T], T]
ValueOrObjectArray = Union[List[ObjectWithIdAndValue[T]], List[T]]


def is_object_with_id_and_value(value_or_object: ValueOrObject[T]) -> bool:
    return isinstance(value_or_object, ObjectWithIdAndValue)


class ScalarWithId(Generic[T]):
    def __init__(self, value_or_object: ValueOrObject[T], id: int = 0):
        if is_object_with_id_and_value(value_or_object):
            self.id = value_or_object.id
            self.value = value_or_object.value
        else:
            self.id = id
            self.value = value_or_object

    def to_json(self) -> Dict[str, Union[int, T]]:
        return {
            "id": self.id,
            "value": self.value
        }
