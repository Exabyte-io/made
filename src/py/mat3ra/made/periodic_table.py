import json
from pathlib import Path
from typing import Dict, Optional

from mat3ra.periodic_table import PERIODIC_TABLE


_PERIODIC_TABLE_DATA: Optional[Dict[str, float]] = None


def _load_periodic_table_data() -> Dict[str, float]:
    global _PERIODIC_TABLE_DATA
    if _PERIODIC_TABLE_DATA is not None:
        return _PERIODIC_TABLE_DATA

    if PERIODIC_TABLE is None:
        raise ImportError(
            "mat3ra.periodic_table is required for periodic table functionality. "
            "Install it with: pip install mat3ra-made[tools] or pip install mat3ra-periodic-table"
        )

    _PERIODIC_TABLE_DATA = {}
    possible_mass_fields = [
        "atomic_mass",
        "atomicMass",
        "mass",
        "atomic_weight",
        "atomicWeight",
        "standard_atomic_weight",
    ]

    for symbol, element_data in PERIODIC_TABLE.items():
        atomic_mass = None
        for field in possible_mass_fields:
            if field in element_data:
                atomic_mass = element_data[field]
                break

        if atomic_mass is not None:
            _PERIODIC_TABLE_DATA[symbol] = float(atomic_mass)
        else:
            available_keys = list(element_data.keys())
            raise ValueError(
                f"Atomic mass not found for element {symbol} in PERIODIC_TABLE. " f"Available keys: {available_keys}"
            )

    return _PERIODIC_TABLE_DATA


def get_atomic_mass_from_element(element: str) -> float:
    """
    Get atomic mass for an element symbol.

    Args:
        element: Element symbol (e.g., "Si", "H", "O")

    Returns:
        Atomic mass in atomic mass units (amu)

    Raises:
        ValueError: If element symbol is not found
        ImportError: If mat3ra.periodic_table is not installed
    """
    data = _load_periodic_table_data()
    element_upper = element.strip().capitalize()
    if element_upper not in data:
        raise ValueError(f"Element symbol '{element}' not found in periodic table")
    return data[element_upper]


def export_periodic_table_json(output_path: Optional[Path] = None) -> Dict[str, float]:
    """
    Export periodic table data as JSON.

    Args:
        output_path: Optional path to save JSON file. If None, returns dict only.

    Returns:
        Dictionary mapping element symbols to atomic masses
    """
    data = _load_periodic_table_data()
    if output_path is not None:
        with open(output_path, "w") as f:
            json.dump(data, f, indent=2)
    return data
