from typing import Dict, Union, List, Tuple, TypedDict
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder, ZSLGenerator
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
import numpy as np

from .interface import InterfaceDataHolder, StrainModes, patch_interface_with_mean_abs_strain
from ..utils import translate_to_bottom_pymatgen_structure
from ..convert import decorator_convert_material_args_kwargs_to_structure


class Settings(TypedDict):
    SUBSTRATE_PARAMETERS: TypedDict[Tuple[int, int, int], int]
    LAYER_PARAMETERS: TypedDict[Tuple[int, int, int], int]
    USE_CONVENTIONAL_CELL: bool
    ZSL_PARAMETERS: TypedDict[float, float, float, float]
    INTERFACE_PARAMETERS: TypedDict[str, Union[str, int, float, tuple]]


@decorator_convert_material_args_kwargs_to_structure
def create_interfaces(
    substrate: Structure,
    layer: Structure,
    settings: Settings,
    sort_by_strain_and_size: bool = True,
    remove_duplicates: bool = True,
) -> InterfaceDataHolder:
    """
    Create all interfaces between the substrate and layer structures using ZSL algorithm provided by pymatgen.

    Args:
        substrate (Structure): The substrate structure.
        layer (Structure): The layer structure.
        settings: The settings for the interface generation.
        sort_by_strain_and_size (bool): Whether to sort the interfaces by strain and size.
        remove_duplicates (bool): Whether to remove duplicate interfaces.
    Returns:
        Dict[str, List[Dict[str, Union[Structure, np.ndarray]]]]: A dictionary of interfaces for each
        termination.
    """
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
            sort_interfaces_by_strain_and_size=sort_by_strain_and_size,
            remove_duplicates=remove_duplicates,
        )
    unique_str = "unique" if remove_duplicates else ""
    print(f"Found {len(interfaces_data.get_interfaces_for_termination(0))} {unique_str} interfaces.")
    return interfaces_data


def normalize_structure(structure: Structure, use_conventional_cell: bool = True):
    """
    Translate atoms to the bottom of the cell (vacuum on top) to allow for the correct consecutive interface generation.
    If use_conventional_cell is passed, conventional cell is used.

    Args:
        structure (Structure): The pymatgen Structure object to normalize.
        use_conventional_cell: Whether to convert to the conventional cell.
    Returns:
        Structure: The normalized pymatgen Structure object.
    """
    if use_conventional_cell:
        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
    structure = translate_to_bottom_pymatgen_structure(structure)
    return structure


@decorator_convert_material_args_kwargs_to_structure
def interface_init_zsl_builder(substrate: Structure, layer: Structure, settings: Settings) -> CoherentInterfaceBuilder:
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
