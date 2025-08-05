def create_point_defect_vacancy(
    material: Union[Material, MaterialWithBuildMetadata],
    coordinate: List[float],
    placement_method: str,
    use_cartesian_coordinates: bool = False,
) -> Material:
    """
    Create a vacancy defect in the given material at the specified coordinate.

    Args:
        material (Material): The host material.
        coordinate (List[float]): The coordinate where the vacancy will be created.
        placement_method (VacancyPlacementMethodEnum): Method to resolve the final coordinate.
        use_cartesian_coordinates (bool): Whether the input coordinate is in Cartesian units.

    Returns:
        Material: A new material with the vacancy defect.
    """
    if placement_method not in [e.value for e in VacancyPlacementMethodEnum]:
        raise ValueError(f"Unsupported placement method: {placement_method}")

    # Convert coordinate to crystal if needed
    if use_cartesian_coordinates:
        coordinate = material.basis.cell.convert_point_to_crystal(coordinate)

    analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
    resolved_coordinate = analyzer.closest_site_coordinate
    config = VacancyDefectConfiguration.from_parameters(crystal=material, coordinate=resolved_coordinate)
    builder = VacancyDefectBuilder()
    return builder.get_material(config)
