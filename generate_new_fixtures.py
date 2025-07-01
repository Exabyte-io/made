#!/usr/bin/env python3
"""
Script to generate new fixture data for nanoribbon tests.
"""

import json
import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src", "py"))

from mat3ra.made.material import Material
from mat3ra.made.tools.build.nanoribbon.helpers import create_nanoribbon

# Test data from the test file
GRAPHENE = {
    "basis": {
        "constraints": [],
        "coordinates": [
            {"id": 0, "value": [0, 0, 0]},
            {"id": 1, "value": [0.333333, 0.666667, 0]},
        ],
        "elements": [
            {"id": 0, "value": "C"},
            {"id": 1, "value": "C"},
        ],
        "labels": [],
        "units": "crystal",
    },
    "isNonPeriodic": False,
    "lattice": {
        "a": 2.467291,
        "alpha": 90,
        "b": 2.467291,
        "beta": 90,
        "c": 20.0,
        "gamma": 120,
        "type": "HEX",
        "units": {"angle": "degree", "length": "angstrom"},
    },
    "name": "Graphene",
}


def generate_new_fixtures():
    """Generate new fixture data for the nanoribbon tests."""
    material = Material.create(GRAPHENE)

    # Generate zigzag nanoribbon
    zigzag_nanoribbon = create_nanoribbon(
        material=material,
        miller_indices_2d=(1, 1),  # zigzag
        width=2,
        length=4,
        vacuum_width=10.0,
        vacuum_length=0.0,
    )

    # Generate armchair nanoribbon
    armchair_nanoribbon = create_nanoribbon(
        material=material,
        miller_indices_2d=(0, 1),  # armchair
        width=2,
        length=4,
        vacuum_width=10.0,
        vacuum_length=0.0,
    )

    # Convert to JSON
    zigzag_json = json.loads(zigzag_nanoribbon.to_json())
    armchair_json = json.loads(armchair_nanoribbon.to_json())

    print("GRAPHENE_ZIGZAG_NANORIBBON = {")
    print(json.dumps(zigzag_json, indent=4))
    print("}")
    print()
    print("GRAPHENE_ARMCHAIR_NANORIBBON = {")
    print(json.dumps(armchair_json, indent=4))
    print("}")


if __name__ == "__main__":
    generate_new_fixtures()
