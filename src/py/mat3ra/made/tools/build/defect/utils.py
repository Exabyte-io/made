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
        # Try direct conversion first
        try:
            return enum_class(value)
        except ValueError:
            # If direct conversion fails, try to find a matching value
            for enum_value in enum_class:
                if enum_value.value == value:
                    return enum_value
            # If no match found, raise the original error with helpful message
            valid_values = [e.value for e in enum_class]
            raise ValueError(f"'{value}' is not a valid {enum_class.__name__}. Valid values are: {valid_values}")
    return value
