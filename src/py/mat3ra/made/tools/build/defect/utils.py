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
        # Get all valid values, handling both direct values and referenced enum values
        valid_values = {}
        for enum_value in enum_class:
            actual_value = get_enum_value(enum_value)
            valid_values[actual_value] = enum_value
        
        # Try to find a matching value
        if value in valid_values:
            return valid_values[value]
            
        # If no match found, raise error with helpful message
        raise ValueError(
            f"'{value}' is not a valid {enum_class.__name__}. "
            f"Valid values are: {list(valid_values.keys())}"
        )
    return value
