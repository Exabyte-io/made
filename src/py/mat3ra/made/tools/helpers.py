"""
Made Tools - Centralized Helper Functions

This module provides convenient access to all helper functions across the Made Tools package.
Functions are organized by category for better discoverability and maintainability.
"""

# Surface Defects
from .build.defect.adatom.helpers import create_adatom_defect, get_adatom_defect_analyzer_cls
from .build.defect.island.helpers import create_island_defect, get_coordinate_condition

# Complex Defects
from .build.defect.pair_defect.helpers import create_pair_defect

# ===== DEFECT FUNCTIONS =====
# Point Defects
from .build.defect.point.helpers import (
    create_point_defect_interstitial,
    create_point_defect_substitution,
    create_point_defect_vacancy,
)
from .build.defect.slab.helpers import create_slab_stack, recreate_slab_with_fractional_layers
from .build.defect.terrace.helpers import create_terrace

# ===== GRAIN BOUNDARY FUNCTIONS =====
from .build.grain_boundary.helpers import create_grain_boundary_linear, create_grain_boundary_planar

# ===== INTERFACE FUNCTIONS =====
# Base Interface
from .build.interface.base.helpers import create_simple_interface_between_slabs

# Commensurate Interfaces
from .build.interface.commensurate.helpers import (
    create_commensurate_interface,
    get_commensurate_strained_configurations,
)

# Twisted Interfaces
from .build.interface.twisted.helpers import create_twisted_interface

# ===== UTILITY FUNCTIONS =====
# Interface Utilities
from .build.interface.utils import get_optimal_film_displacement, get_slab, remove_duplicate_interfaces

# ZSL Interfaces
from .build.interface.zsl.helpers import create_zsl_interface, create_zsl_interface_between_slabs

# Lattice Utilities
from .build.lattice_lines.helpers import create_lattice_lines_config_and_material
from .build.monolayer.helpers import create_monolayer

# ===== NANOSTRUCTURE FUNCTIONS =====
from .build.nanoparticle.helpers import create_nanoparticle_from_material
from .build.nanoribbon.helpers import create_nanoribbon
from .build.nanotape.helpers import create_nanotape

# ===== PASSIVATION FUNCTIONS =====
from .build.passivation.helpers import get_unique_coordination_numbers, passivate_surface
from .build.perturbation.helpers import create_perturbation

# ===== CORE BUILD FUNCTIONS =====
from .build.slab.helpers import create_atomic_layers, create_slab, get_slab_terminations
from .build.supercell.helpers import create_supercell

# ===== EXPORT ALL =====
__all__ = [
    # Core Build Functions
    "create_slab",
    "create_atomic_layers",
    "get_slab_terminations",
    "create_supercell",
    "create_monolayer",
    "create_perturbation",
    # Nanostructure Functions
    "create_nanoparticle_from_material",
    "create_nanoribbon",
    "create_nanotape",
    # Point Defect Functions
    "create_point_defect_vacancy",
    "create_point_defect_substitution",
    "create_point_defect_interstitial",
    # Surface Defect Functions
    "create_adatom_defect",
    "get_adatom_defect_analyzer_cls",
    "create_island_defect",
    "get_coordinate_condition",
    "create_terrace",
    # Complex Defect Functions
    "create_pair_defect",
    "create_slab_stack",
    "recreate_slab_with_fractional_layers",
    # Interface Functions
    "create_simple_interface_between_slabs",
    "create_zsl_interface",
    "create_zsl_interface_between_slabs",
    "create_commensurate_interface",
    "get_commensurate_strained_configurations",
    "create_twisted_interface",
    # Grain Boundary Functions
    "create_grain_boundary_planar",
    "create_grain_boundary_linear",
    # Passivation Functions
    "passivate_surface",
    "get_unique_coordination_numbers",
    # Utility Functions
    "get_optimal_film_displacement",
    "get_slab",
    "remove_duplicate_interfaces",
    "create_lattice_lines_config_and_material",
]
