from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder, ZSLGenerator
from pymatgen.core.structure import Structure
import numpy as np

from .modify import translate_to_bottom

strain_modes_map = {
    "strain": "strain",
    "von_mises_strain": "von_mises_strain",
    "mean_abs_strain": "mean_abs_strain",
}


# TODO: add decorator to convert ESSE Material to Pymatgen
def create_interfaces(substrate: Structure, layer: Structure, settings):
    translate_to_bottom(substrate)
    translate_to_bottom(layer)
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

    # Find terminations
    cib._find_terminations()
    terminations = cib.terminations

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
            # Wrap atoms to unit cell
            interface.make_supercell((1, 1, 1), to_unit_cell=True)
            mean_abs_strain = round(
                np.mean(np.abs(interface.interface_properties["strain"])) / settings["INTERFACE_PARAMETERS"]["STRAIN_TOLERANCE"]
            ) * settings["INTERFACE_PARAMETERS"]["STRAIN_TOLERANCE"]
            interface_struct = {
                "interface": interface,
                "strain": interface.interface_properties["strain"],
                "von_mises_strain": interface.interface_properties["strain"].von_mises_strain,
                "mean_abs_strain": mean_abs_strain
            }
            interfaces[termination].append(interface_struct)
    return interfaces, terminations


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
