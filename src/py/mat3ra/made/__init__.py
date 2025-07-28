import os
import warnings
from pathlib import Path


def _is_development_mode():
    """
    Determine if we're running in development mode.
    Returns True if we're in development, False if installed as a package.
    """
    # Check if we're in a git repository (development mode indicator)
    current_file = Path(__file__).resolve()
    
    # Look for .git directory or setup.py/pyproject.toml in parent directories
    for parent in current_file.parents:
        if (parent / ".git").exists() or (parent / "pyproject.toml").exists():
            return True
        # Stop at root directory
        if parent == parent.parent:
            break
    
    # Check if MADE_DEV_MODE environment variable is set
    return os.environ.get("MADE_DEV_MODE", "").lower() in ("1", "true", "yes")


def configure_warnings(show_pydantic_warnings=None):
    """
    Configure warnings for the mat3ra-made package.
    
    Args:
        show_pydantic_warnings (bool, optional): 
            If True, show pydantic warnings.
            If False, suppress pydantic warnings.
            If None, auto-detect based on environment (show in dev, hide in production).
    """
    if show_pydantic_warnings is None:
        # Auto-detect: show warnings in development, hide in production
        show_pydantic_warnings = _is_development_mode()
    
    if not show_pydantic_warnings:
        # Suppress pydantic serialization warnings
        warnings.filterwarnings(
            "ignore", 
            category=UserWarning, 
            module="pydantic.*",
            message=".*Pydantic serializer warnings.*"
        )
        warnings.filterwarnings(
            "ignore", 
            category=UserWarning, 
            module="pydantic.*",
            message=".*serialized value may not be as expected.*"
        )


def suppress_pydantic_warnings():
    """
    Convenience function to suppress pydantic warnings.
    Useful for notebook environments where you want to explicitly control warnings.
    """
    configure_warnings(show_pydantic_warnings=False)


def enable_pydantic_warnings():
    """
    Convenience function to enable pydantic warnings.
    Useful for development or debugging.
    """
    configure_warnings(show_pydantic_warnings=True)


# Auto-configure warnings on import based on environment
configure_warnings()

# Public API
__all__ = [
    "configure_warnings",
    "suppress_pydantic_warnings", 
    "enable_pydantic_warnings"
]
