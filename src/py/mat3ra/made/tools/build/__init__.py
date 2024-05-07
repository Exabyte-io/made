from ...material import Material
from .interface import InterfaceDataHolder, CoherentInterfaceBuilder, TerminationType
from .interface import InterfaceSettings as Settings
from .interface import interface_init_zsl_builder, interface_patch_with_mean_abs_strain
from ..convert import decorator_convert_material_args_kwargs_to_structure
from ..modify import translate_to_bottom, wrap_to_unit_cell


@decorator_convert_material_args_kwargs_to_structure
def create_interfaces(
    substrate: Material = None,
    layer: Material = None,
    settings: Settings = None,
    sort_by_strain_and_size: bool = True,
    remove_duplicates: bool = True,
    is_logging_enabled: bool = True,
    interface_builder: CoherentInterfaceBuilder = None,
    termination: TerminationType = None,
) -> InterfaceDataHolder:
    """
    Create all interfaces between the substrate and layer structures using ZSL algorithm provided by pymatgen.

    Args:
        substrate (Material): The substrate structure.
        layer (Material): The layer structure.
        settings: The settings for the interface generation.
        sort_by_strain_and_size (bool): Whether to sort the interfaces by strain and size.
        remove_duplicates (bool): Whether to remove duplicate interfaces.
        is_logging_enabled (bool): Whether to enable debug print.
    Returns:
        InterfaceDataHolder.
    """
    if is_logging_enabled:
        print("Creating interfaces...")

    builder = interface_builder or init_interface_builder(substrate, layer, settings)
    interfaces_data = InterfaceDataHolder()

    if termination is not None:
        builder.terminations = [termination]

    for termination in builder.terminations:
        all_interfaces_for_termination = builder.get_interfaces(
            termination,
            gap=settings["INTERFACE_PARAMETERS"]["DISTANCE_Z"],
            film_thickness=settings["LAYER_PARAMETERS"]["THICKNESS"],
            substrate_thickness=settings["SUBSTRATE_PARAMETERS"]["THICKNESS"],
            in_layers=True,
        )

        all_interfaces_for_termination_patched_wrapped = list(
            map(
                lambda i: wrap_to_unit_cell(interface_patch_with_mean_abs_strain(i)),
                all_interfaces_for_termination,
            )
        )

        interfaces_data.add_data_entries(
            all_interfaces_for_termination_patched_wrapped,
            sort_interfaces_by_strain_and_size=sort_by_strain_and_size,
            remove_duplicates=remove_duplicates,
        )

    if is_logging_enabled:
        unique_str = "unique" if remove_duplicates else ""
        print(f"Found {len(interfaces_data.get_interfaces_for_termination(0))} {unique_str} interfaces.")

    return interfaces_data


@decorator_convert_material_args_kwargs_to_structure
def init_interface_builder(
    substrate: Material,
    layer: Material,
    settings: Settings,
) -> CoherentInterfaceBuilder:
    substrate = translate_to_bottom(substrate, settings["USE_CONVENTIONAL_CELL"])
    layer = translate_to_bottom(layer, settings["USE_CONVENTIONAL_CELL"])
    return interface_init_zsl_builder(substrate, layer, settings)
