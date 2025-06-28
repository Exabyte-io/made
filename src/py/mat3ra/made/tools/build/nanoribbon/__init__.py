# Convenience imports for the nanoribbon package
from .configuration import (
    CrystalLatticeLinesConfiguration,
    CrystalLatticeLinesUniqueRepeatedConfiguration,
    NanoribbonConfiguration,
    get_miller_indices_from_edge_type,
)

# Import only specific builders that don't conflict
from .builders import (
    CrystalLatticeLinesBuilder,
    CrystalLatticeLinesRepeatedBuilder,
)

# Keep the enums for backward compatibility
from .enums import EdgeTypes

# Import helper functions explicitly to avoid conflicts
from .helpers import (
    create_nanoribbon as _create_nanoribbon_function,
    create_nanoribbon_from_edge_type,
    get_nanoribbon_terminations,
)

# Explicitly assign to ensure correct function is used
create_nanoribbon = _create_nanoribbon_function

__all__ = [
    # Configurations
    "CrystalLatticeLinesConfiguration",
    "CrystalLatticeLinesUniqueRepeatedConfiguration", 
    "NanoribbonConfiguration",
    # Builders (NanoribbonBuilder removed to avoid method conflicts)
    "CrystalLatticeLinesBuilder",
    "CrystalLatticeLinesRepeatedBuilder",
    # Helper functions
    "create_nanoribbon",
    "create_nanoribbon_from_edge_type",
    "get_nanoribbon_terminations",
    "get_miller_indices_from_edge_type",
    # Legacy
    "EdgeTypes",
]

# Note: NanoribbonBuilder can still be imported directly as:
# from mat3ra.made.tools.build.nanoribbon.builders import NanoribbonBuilder 