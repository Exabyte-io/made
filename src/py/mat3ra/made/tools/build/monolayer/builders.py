from typing import Optional

from mat3ra.made.material import Material
from .configurations import MonolayerConfiguration
from .. import BaseSingleBuilder
from ..slab.helpers import create_slab
from ...modify import translate_to_z_level, filter_by_box, add_vacuum


class MonolayerBuilder(BaseSingleBuilder):
    """
    Builder for creating monolayer structures from crystal materials.
    
    The builder creates different monolayer structures based on the crystal type:
    - HEX: Creates a slab with Miller indices (0,0,1), thickness=1, translates to center, then filters half
    - FCC/CUB: Creates a slab with Miller indices (1,1,1), thickness=1 using primitive cell, then applies
      specific filtering and centering operations
    """
    
    _ConfigurationType = MonolayerConfiguration
    _GeneratedItemType = Material
    
    def _generate(self, configuration: MonolayerConfiguration) -> Material:
        crystal = configuration.crystal
        lattice_type_str = crystal.lattice.type.value if hasattr(crystal.lattice.type, 'value') else str(crystal.lattice.type)
        vacuum = configuration.vacuum
        
        if lattice_type_str == "HEX":
            return self._create_hex_monolayer(crystal, vacuum)
        elif lattice_type_str in ["FCC", "CUB"]:
            return self._create_fcc_cub_monolayer(crystal, vacuum)
        else:
            return self._create_generic_monolayer(crystal, vacuum)

    
    def _create_hex_monolayer(
        self, 
        crystal: Material, 
        vacuum: float = 0.0
    ) -> Material:
        miller_indices = (0, 0, 1)
        
        slab = create_slab(
            crystal=crystal,
            miller_indices=miller_indices,
            number_of_layers=1,
            vacuum=vacuum,
            use_conventional_cell=True,
        )
        
        centered_slab = translate_to_z_level(slab, z_level="center")
        
        half_filtered_slab = filter_by_box(
            centered_slab, 
            min_coordinate=[0, 0, 0.5], 
            max_coordinate=[1, 1, 1], 
            use_cartesian_coordinates=False
        )

            
        return half_filtered_slab
    
    def _create_fcc_cub_monolayer(
        self, 
        crystal: Material, 
        vacuum: float
    ) -> Material:
        miller_indices = (1, 1, 1)
        
        slab = create_slab(
            crystal=crystal,
            miller_indices=miller_indices,
            number_of_layers=1,
            vacuum=vacuum,
            use_conventional_cell=False,
        )
        
        return slab
    
    def _create_generic_monolayer(
        self, 
        crystal: Material, 
        vacuum: float
    ) -> Material:
        return self._create_fcc_cub_monolayer(crystal, vacuum)
    
    def _update_material_name(self, material: Material, configuration: MonolayerConfiguration) -> Material:
        crystal = configuration.crystal
        lattice_type_str = crystal.lattice.type.value if hasattr(crystal.lattice.type, 'value') else str(crystal.lattice.type)
        
        if lattice_type_str == "HEX":
            miller_str = "001"
        elif lattice_type_str in ["FCC", "CUB"] or self._is_cubic_like(crystal):
            miller_str = "111"
        else:
            miller_str = "111"
            
        original_name = crystal.name or "Crystal"
        material.name = f"{original_name} - Monolayer ({miller_str})"
        return material 
