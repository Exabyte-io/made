from typing import Tuple, Dict, Union
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder, ZSLGenerator
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
import numpy as np

from .utils import translate_to_bottom_pymatgen_structure
from .convert import convert_material_args_kwargs_to_structure

strain_modes_map = {
    "strain": "strain",
    "von_mises_strain": "von_mises_strain",
    "mean_abs_strain": "mean_abs_strain",
}


class InterfaceSettings:
    SUBSTRATE_PARAMETERS: Dict[str, Union[Tuple[int, int, int], int]] = {
        "MILLER_INDICES": (1, 1, 1),  # the miller indices of the interfacial plane
        "THICKNESS": 3,  # in layers
    }
    LAYER_PARAMETERS: Dict[str, Union[Tuple[int, int, int], int]] = {
        "MILLER_INDICES": (0, 0, 1),  # the miller indices of the interfacial plane
        "THICKNESS": 1,  # in layers
    }
    USE_CONVENTIONAL_CELL: bool = True
    ZSL_PARAMETERS: Dict[str, float] = {
        "MAX_AREA": 50,  # The area to consider in Angstrom^2
        "MAX_AREA_TOL": 0.09,  # The area within this tolerance is considered equal
        "MAX_LENGTH_TOL": 0.03,  # supercell lattice vectors lengths within this tolerance are considered equal
        "MAX_ANGLE_TOL": 0.01,  # supercell lattice angles within this tolerance are considered equal
        "STRAIN_TOL": 10e-6,  # strains within this tolerance are considered equal
    }
    INTERFACE_PARAMETERS: Dict[str, float] = {
        "DISTANCE_Z": 3.0,  # in Angstroms
        "MAX_AREA": 50,  # in Angstroms^2
    }


@convert_material_args_kwargs_to_structure
def create_interfaces(substrate: Structure, layer: Structure, settings: InterfaceSettings):
    """
    Create all interfaces between the substrate and layer structures using ZSL algorithm provided by pymatgen.
    Args:
        substrate (Structure): The substrate structure.
        layer (Structure): The layer structure.
        settings (InterfaceSettings): The settings for the interface generation.

    Returns:
        Tuple[Dict[str, List[Dict[str, Union[Structure, np.ndarray]]]], List[str]]: A dictionary of interfaces for each
        termination and a list of terminations.
    """
    substrate = normalize_structure(substrate, settings["USE_CONVENTIONAL_CELL"])
    layer = normalize_structure(layer, settings["USE_CONVENTIONAL_CELL"])

    print("Creating interfaces...")
    zsl: ZSLGenerator = ZSLGenerator(
        max_area_ratio_tol=settings["ZSL_PARAMETERS"]["MAX_AREA_TOL"],
        max_area=settings["ZSL_PARAMETERS"]["MAX_AREA"],
        max_length_tol=settings["ZSL_PARAMETERS"]["MAX_LENGTH_TOL"],
        max_angle_tol=settings["ZSL_PARAMETERS"]["MAX_ANGLE_TOL"],
    )

    cib: CoherentInterfaceBuilder = CoherentInterfaceBuilder(
        substrate_structure=substrate,
        film_structure=layer,
        substrate_miller=settings["SUBSTRATE_PARAMETERS"]["MILLER_INDICES"],
        film_miller=settings["LAYER_PARAMETERS"]["MILLER_INDICES"],
        zslgen=zsl,
    )

    terminations = cib.terminations
    print(terminations)

    # Function: create interfaces for terminations
    # Create interfaces for each termination
    interfaces = {}
    for termination in terminations:
        interfaces[termination] = []
        all_interfaces_for_termination = cib.get_interfaces(
            termination,
            gap=settings["INTERFACE_PARAMETERS"]["DISTANCE_Z"],
            film_thickness=settings["LAYER_PARAMETERS"]["THICKNESS"],
            substrate_thickness=settings["SUBSTRATE_PARAMETERS"]["THICKNESS"],
            in_layers=True,
        )

        for interface in all_interfaces_for_termination:
            # Function: add data to interface
            # Wrap atoms to unit cell
            interface.make_supercell((1, 1, 1), to_unit_cell=True)
            mean_abs_strain = (
                round(
                    np.mean(np.abs(interface.interface_properties["strain"])) / settings["ZSL_PARAMETERS"]["STRAIN_TOL"]
                )
                * settings["ZSL_PARAMETERS"]["STRAIN_TOL"]
            )
            interface_struct = {
                "interface": interface,
                "strain": interface.interface_properties["strain"],
                "von_mises_strain": interface.interface_properties["von_mises_strain"],
                "mean_abs_strain": mean_abs_strain,
            }
            interfaces[termination].append(interface_struct)
    return interfaces, terminations


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
    structure = translate_to_bottom_pymatgen_structure(structure)
    if conventional_cell:
        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
    return structure


def get_terminations(pymatgen_materials, settings):
    zsl: ZSLGenerator = ZSLGenerator(
        max_area_ratio_tol=settings["ZSL_PARAMETERS"]["MAX_AREA_TOL"],
        max_area=settings["ZSL_PARAMETERS"]["MAX_AREA"],
        max_length_tol=settings["ZSL_PARAMETERS"]["MAX_LENGTH_TOL"],
        max_angle_tol=settings["ZSL_PARAMETERS"]["MAX_ANGLE_TOL"],
    )

    cib = CoherentInterfaceBuilder(
        substrate_structure=pymatgen_materials[settings["SUBSTRATE_PARAMETERS"]["MATERIAL_INDEX"]],
        film_structure=pymatgen_materials[settings["LAYER_PARAMETERS"]["MATERIAL_INDEX"]],
        substrate_miller=settings["SUBSTRATE_PARAMETERS"]["MILLER_INDICES"],
        film_miller=settings["LAYER_PARAMETERS"]["MILLER_INDICES"],
        zslgen=zsl,
    )

    # Find terminations
    cib._find_terminations()
    terminations = cib.terminations
    return terminations


def sort_interfaces(interfaces, terminations, strain_mode=strain_modes_map["mean_abs_strain"]):
    sorted_interfaces = {}
    for termination in terminations:
        sorted_interfaces[termination] = sorted(
            interfaces[termination], key=lambda x: (x[strain_mode], x["interface"].num_sites)
        )
    return sorted_interfaces
