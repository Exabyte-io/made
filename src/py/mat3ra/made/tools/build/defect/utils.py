from enum import Enum
from typing import TypeVar, Type, Union

EnumType = TypeVar("EnumType", bound=Enum)


def ensure_enum(value: Union[str, EnumType], enum_class: Type[EnumType]) -> EnumType:
    """
    Ensures a value is converted to the specified enum type if it's a string.

    Args:
        value: String or enum value
        enum_class: The enum class to convert to

    Returns:
        EnumType: The enum value
    """
    if isinstance(value, str):
        return enum_class(value)
    return value
