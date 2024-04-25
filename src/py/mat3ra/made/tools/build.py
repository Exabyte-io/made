from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder, ZSLGenerator
from pymatgen.core.interface import Interface

yield {
    "interface": Interface.from_slabs(
        substrate_slab=sub_sl_slab,
        film_slab=film_sl_slab,
        gap=gap,
        vacuum_over_film=vacuum_over_film,
        interface_properties=interface_properties,
        center_slab=False,  # False -- positions interface at the most bottom of the cell, solving the issue of second iteration not working properly
    ),
    "strain": strain,
    "von_mises_strain": strain.von_mises_strain,
    "mean_abs_strain": round(np.mean(np.abs(strain)) / self.strain_tol) * self.strain_tol,
    "film_sl_vectors": match.film_sl_vectors,
    "substrate_sl_vectors": match.substrate_sl_vectors,
    "film_transform": super_film_transform,
    "substrate_transform": super_sub_transform,
}


def create_interfaces(pymatgen_materials, settings):
    print("Creating interfaces...")
    zsl = ZSLGenerator(
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
        strain_tol=settings["ZSL_PARAMETERS"]["STRAIN_TOL"],
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
            interface["interface"].make_supercell((1, 1, 1), to_unit_cell=True)
            interfaces[termination].append(interface)
    return interfaces, terminations


def sort_interfaces(interfaces, terminations, strain_mode="von_mises_strain"):
    sorted_interfaces = {}
    for termination in terminations:
        sorted_interfaces[termination] = sorted(
            interfaces[termination], key=lambda x: (x[strain_mode], x["interface"].num_sites)
        )
    return sorted_interfaces
