def get_passivant_coordinate_crystal(
    material: Material, atom_coordinate: np.ndarray, bond_direction: np.ndarray, bond_length: float = 2.0
):
    bond_normal_cartesian = np.array(material.basis.cell.convert_point_to_cartesian(bond_direction.tolist()))
    bond_normal_cartesian /= np.linalg.norm(bond_normal_cartesian)
    bond_vector_cartesian = bond_normal_cartesian * bond_length
    bond_vector_crystal = np.array(material.basis.cell.convert_point_to_crystal(bond_vector_cartesian.tolist()))
    passivant_coordinate_crystal = atom_coordinate + bond_vector_crystal
    return passivant_coordinate_crystal


def passivate_material(
    slab: Material, passivant: str, bond_length: float = 2.0, coordination_threshold: Optional[int] = None
):
    nudge_value = 0.01
    supercell_scaling_factor = [3, 3, 3]
    min_coordinate = [1 / 3, 1 / 3, 1 / 3]
    max_coordinate = [2 / 3, 2 / 3, 2 / 3]
    adjusted_min_coordinate = (np.array(min_coordinate) - nudge_value).tolist()
    adjusted_max_coordinate = (np.array(max_coordinate) + nudge_value).tolist()
    slab = translate_to_z_level(slab, "center")
    new_basis = slab.basis.copy()
    slab_supercell = create_supercell(slab, scaling_factor=supercell_scaling_factor)

    central_cell_ids = [
        id
        for id, coordinate in zip(slab_supercell.basis.coordinates.ids, slab_supercell.basis.coordinates.values)
        if is_coordinate_in_box(coordinate, adjusted_min_coordinate, adjusted_max_coordinate)
    ]
    undercoordinated_atom_indices, atom_neighbors_info = get_undercoordinated_atoms(
        slab_supercell, central_cell_ids, coordination_threshold
    )
    for index in undercoordinated_atom_indices:
        atom_coordinate_crystal = np.array(slab_supercell.basis.coordinates.values[index])
        neighbors_average_coordinate_crystal = atom_neighbors_info[index][2]
        bond_normal_crystal = np.array(slab_supercell.basis.coordinates.values[index]) - np.array(
            neighbors_average_coordinate_crystal
        )
        passivant_coordinate_crystal = get_passivant_coordinate_crystal(
            slab_supercell, atom_coordinate_crystal, bond_normal_crystal, bond_length
        ).tolist()
        passivant_coordinate_crystal_original_cell = transform_coordinate_to_supercell(
            passivant_coordinate_crystal,
            scaling_factor=supercell_scaling_factor,
            translation_vector=min_coordinate,
            reverse=True,
        )
        new_basis.add_atom(passivant, passivant_coordinate_crystal_original_cell)

    slab.basis = new_basis
    return slab


def add_passivant_atoms_to_material(
    material: Material,
    axis: Literal["x", "y", "z"],
    surface: Optional[SURFACE_TYPES] = SURFACE_TYPES.BOTH,
    tolerance: float = 0.01,
    passivant: Optional[str] = "H",
    default_bond_length: float = 1.0,
) -> Material:
    """
    Add passivant atoms to the specified surface or edge of the material.

    Args:
        material (Material): The material object to add passivant atoms to.
        axis (AXIS_TYPES): The axis along which the surface or edge lies ("x", "y", or "z").
        surface (SURFACE_TYPES): The surface to add passivant atoms to ("top", "bottom", or "both").
        tolerance (float): The tolerance for selecting surface atoms.
        passivant (str): The chemical symbol of the passivating atom (e.g., 'H').
        default_bond_length (float): The default bond length to use if the pair is not found in BOND_LENGTHS_MAP.

    Returns:
        Material: The material object with passivation atoms added.
    """
    axis_idx = {"x": 0, "y": 1, "z": 2}[axis]
    basis = material.basis.copy()
    if axis == "z":
        if surface in [SURFACE_TYPES.TOP, SURFACE_TYPES.BOTH]:
            z_extremum = get_atomic_coordinates_extremum(material, "max", axis)
            add_passivant_atoms_to_basis(
                material, basis, axis_idx, z_extremum, tolerance, passivant, default_bond_length, positive=True
            )
        if surface in [SURFACE_TYPES.BOTTOM, SURFACE_TYPES.BOTH]:
            z_extremum = get_atomic_coordinates_extremum(material, "min", axis)
            add_passivant_atoms_to_basis(
                material, basis, axis_idx, z_extremum, tolerance, passivant, default_bond_length, positive=False
            )
    else:
        x_or_y_extremum_max = get_atomic_coordinates_extremum(material, "max", axis)
        add_passivant_atoms_to_basis(
            material, basis, axis_idx, x_or_y_extremum_max, tolerance, passivant, default_bond_length, positive=True
        )

        x_or_y_extremum_min = get_atomic_coordinates_extremum(material, "min", axis)
        add_passivant_atoms_to_basis(
            material, basis, axis_idx, x_or_y_extremum_min, tolerance, passivant, default_bond_length, positive=False
        )
    material.basis = basis
    return material


def add_passivant_atoms_to_basis(
    material, basis, axis_idx, extremum, tolerance, passivant, default_bond_length, positive=True
):
    """
    Helper function to add passivant atoms to the specified extremum along a given axis.
    """
    min_box = [0, 0, 0]
    max_box = [1, 1, 1]
    min_box[axis_idx] = extremum - tolerance
    max_box[axis_idx] = extremum + tolerance

    surface_atoms = filter_by_box(material, min_box, max_box)
    surface_indices = surface_atoms.basis.coordinates.ids

    for idx in surface_indices:
        atom_coordinate = material.basis.coordinates.values[idx]
        element = material.basis.elements.values[idx]
        bond_length = BOND_LENGTHS_MAP.get((element, passivant), default_bond_length)

        bond_length_vector = [0, 0, 0]
        bond_length_vector[axis_idx] = bond_length if positive else -bond_length
        bond_length_crystal = material.basis.cell.convert_point_to_crystal(bond_length_vector)

        passivant_coordinate = atom_coordinate + bond_length_crystal
        basis.add_atom(passivant, passivant_coordinate)


def passivate_surface(
    material: Material,
    passivant: str = "H",
    default_bond_length: float = 1.0,
    surface: SURFACE_TYPES = SURFACE_TYPES.BOTH,
) -> Material:
    """
    Passivates the top and/or bottom surfaces of a material by adding atoms along the Z-axis,
    with bond lengths determined by the element and passivant.

    Args:
        material (Material): The material to passivate.
        passivant (str): The chemical symbol of the passivating atom (e.g., 'H').
        default_bond_length (float): The default bond length to use if the pair is not found in BOND_LENGTHS_MAP.
        surface (SURFACE_TYPES): The surface to passivate ("top", "bottom", or "both").

    Returns:
        Material: The passivated material.
    """
    material = translate_to_z_level(material, "center")

    return add_passivant_atoms_to_material(
        material=material,
        axis="z",
        surface=surface,
        passivant=passivant,
        default_bond_length=default_bond_length,
    )


def passivate_edges(
    material: Material,
    passivant: str = "H",
    default_bond_length: float = 1.0,
    axis: Literal["x", "y"] = "y",
) -> Material:
    """
    Passivates the edges of a 2D material by adding atoms along the X or Y axis,
    with bond lengths determined by the element and passivant.

    Args:
        material (Material): The material to passivate.
        passivant (str): The chemical symbol of the passivating atom (e.g., 'H').
        default_bond_length (float): The default bond length to use if the pair is not found in BOND_LENGTHS_MAP.
        axis (AXIS_TYPES): The axis along which the edges lie ("x" or "y").

    Returns:
        Material: The passivated material.
    """
    material = translate_to_z_level(material, "center")

    return add_passivant_atoms_to_material(
        material=material, axis=axis, passivant=passivant, default_bond_length=default_bond_length
    )
