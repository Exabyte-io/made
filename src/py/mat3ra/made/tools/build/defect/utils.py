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
        # Get all valid values, handling both direct values and referenced enum values
        valid_values = {}
        for enum_value in enum_class:
            # Handle both direct string values and referenced enum values
            if isinstance(enum_value.value, str):
                valid_values[enum_value.value] = enum_value
            elif isinstance(enum_value.value, Enum):
                valid_values[enum_value.value.value] = enum_value

        # Try to find a matching value
        if value in valid_values:
            return valid_values[value]

        # If no match found, raise error with helpful message
        raise ValueError(
            f"'{value}' is not a valid {enum_class.__name__}. "
            f"Valid values are: {list(valid_values.keys())}"
        )
    return value
