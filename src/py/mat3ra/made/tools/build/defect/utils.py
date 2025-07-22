from enum import Enum
from typing import TypeVar, Type, Union

EnumType = TypeVar("EnumType", bound=Enum)

def get_enum_value(enum_member: Enum) -> str:
    """
    Gets the actual string value from an enum member, handling nested enum references.
    
    Args:
        enum_member: The enum member to get the value from
        
    Returns:
        str: The actual string value
    """
    value = enum_member.value
    while isinstance(value, Enum):
        value = value.value
    return value

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
        # Try to match by the actual string value
        for enum_value in enum_class:
            # Handle direct string value
            if enum_value.value == value:
                return enum_value
            # Handle referenced enum value
            if isinstance(enum_value.value, Enum) and enum_value.value.value == value:
                return enum_value
        # Try to match by the referenced enum's name (e.g., "CLOSEST_SITE")
        for enum_value in enum_class:
            if isinstance(enum_value.value, Enum) and enum_value.value.name.lower() == value.lower():
                return enum_value
            if enum_value.name.lower() == value.lower():
                return enum_value
        # If no match found, raise error with helpful message
        valid_values = []
        for enum_value in enum_class:
            if isinstance(enum_value.value, Enum):
                valid_values.append(enum_value.value.value)
            else:
                valid_values.append(enum_value.value)
        raise ValueError(
            f"'{value}' is not a valid {enum_class.__name__}. "
            f"Valid values are: {valid_values}"
        )
    return value
