from typing import Dict, Union, List
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder, ZSLGenerator
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
import numpy as np

from .interface import InterfaceDataHolder, StrainModes, patch_interface_with_mean_abs_strain
from ..utils import translate_to_bottom_pymatgen_structure
from ..convert import convert_material_args_kwargs_to_structure


@convert_material_args_kwargs_to_structure
def create_interfaces(substrate: Structure, layer: Structure, settings, **kwargs):
    """
    Create all interfaces between the substrate and layer structures using ZSL algorithm provided by pymatgen.
    Args:
        substrate (Structure): The substrate structure.
        layer (Structure): The layer structure.
        settings: The settings for the interface generation.

    Keyword Args:
        sort_by_strain_and_size (bool): Whether to sort the interfaces by strain and size.
        remove_duplicates (bool): Whether to remove duplicate interfaces.
    Returns:
        Dict[str, List[Dict[str, Union[Structure, np.ndarray]]]]: A dictionary of interfaces for each
        termination.
    """

    sort_by_strain_and_size = kwargs.get("sort_by_strain_and_size", True)
    remove_duplicates = kwargs.get("remove_duplicates", True)

    substrate = normalize_structure(substrate, settings["USE_CONVENTIONAL_CELL"])
    layer = normalize_structure(layer, settings["USE_CONVENTIONAL_CELL"])

    print("Creating interfaces...")
    builder = interface_init_zsl_builder(substrate, layer, settings)

    interfaces_data = InterfaceDataHolder()
    for termination in builder.terminations:
        all_interfaces_for_termination = builder.get_interfaces(
            termination,
            gap=settings["INTERFACE_PARAMETERS"]["DISTANCE_Z"],
            film_thickness=settings["LAYER_PARAMETERS"]["THICKNESS"],
            substrate_thickness=settings["SUBSTRATE_PARAMETERS"]["THICKNESS"],
            in_layers=True,
        )

        all_interfaces_for_termination_patched = list(
            map(patch_interface_with_mean_abs_strain, all_interfaces_for_termination)
        )

        interfaces = []
        for interface in all_interfaces_for_termination_patched:
            # Wrap atoms to unit cell
            interface.make_supercell((1, 1, 1), to_unit_cell=True)
            interfaces.append(interface)
        interfaces_data.add_data_entries(
            interfaces,
            sort_interfaces_for_all_terminations_by_strain_and_size=sort_by_strain_and_size,
            remove_duplicates=remove_duplicates,
        )
    unique_str = "unique" if remove_duplicates else ""
    print(f"Found {len(interfaces_data.get_interfaces_for_termination(0))} {unique_str} interfaces.")
    return interfaces_data


def normalize_structure(structure: Structure, conventional_cell: bool = True):
    """
    Normalize the structure to the conventional cell.
    And translate the structure to the bottom of the cell to allow for the correct consecutive interface generation.
    Args:
        structure (Structure): The pymatgen Structure object to normalize.
        conventional_cell: Whether to convert to the conventional cell.
    Returns:
        Structure: The normalized pymatgen Structure object.
    """
    if conventional_cell:
        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
    structure = translate_to_bottom_pymatgen_structure(structure)
    return structure


@convert_material_args_kwargs_to_structure
def interface_init_zsl_builder(substrate: Structure, layer: Structure, settings) -> CoherentInterfaceBuilder:
    generator: ZSLGenerator = ZSLGenerator(
        max_area_ratio_tol=settings["ZSL_PARAMETERS"]["MAX_AREA_TOL"],
        max_area=settings["ZSL_PARAMETERS"]["MAX_AREA"],
        max_length_tol=settings["ZSL_PARAMETERS"]["MAX_LENGTH_TOL"],
        max_angle_tol=settings["ZSL_PARAMETERS"]["MAX_ANGLE_TOL"],
    )

    builder = CoherentInterfaceBuilder(
        substrate_structure=substrate,
        film_structure=layer,
        substrate_miller=settings["SUBSTRATE_PARAMETERS"]["MILLER_INDICES"],
        film_miller=settings["LAYER_PARAMETERS"]["MILLER_INDICES"],
        zslgen=generator,
    )

    return builder


def sort_interfaces(interfaces, terminations, strain_mode=StrainModes.mean_abs_strain):
    sorted_interfaces = {}
    for termination in terminations:
        sorted_interfaces[termination] = sorted(
            interfaces[termination], key=lambda x: (x[strain_mode], x["interface"].num_sites)
        )
    return sorted_interfaces
