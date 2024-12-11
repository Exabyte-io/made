from typing import List, Callable, Optional, Union

import numpy as np
from pydantic import BaseModel
from mat3ra.made.material import Material


from mat3ra.made.utils import get_center_of_coordinates
from ...third_party import (
    PymatgenStructure,
    PymatgenPeriodicSite,
    PymatgenVacancy,
    PymatgenSubstitution,
    PymatgenInterstitial,
)

from ...modify import (
    add_vacuum,
    filter_by_ids,
    filter_by_box,
    filter_by_condition_on_coordinates,
    translate_to_z_level,
    rotate,
)
from ...build import BaseBuilder
from ...convert import to_pymatgen
from ...analyze.other import (
    get_atomic_coordinates_extremum,
    get_closest_site_id_from_coordinate,
    get_closest_site_id_from_coordinate_and_element,
    get_local_extremum_atom_index,
)
from ...analyze.coordination import get_voronoi_nearest_neighbors_atom_indices
from ...utils import transform_coordinate_to_supercell, coordinate as CoordinateCondition
from ..utils import merge_materials
from ..slab import SlabConfiguration, create_slab, Termination
from ..supercell import create_supercell
from ..mixins import ConvertGeneratedItemsPymatgenStructureMixin
from .configuration import (
    PointDefectConfiguration,
    AdatomSlabPointDefectConfiguration,
    IslandSlabDefectConfiguration,
    TerraceSlabDefectConfiguration,
    PointDefectPairConfiguration,
)
from .factories import DefectBuilderFactory


class PointDefectBuilderParameters(BaseModel):
    center_defect: bool = False


class DefectBuilder(BaseBuilder):
    def create_isolated_defect(self, defect_configuration: PointDefectConfiguration) -> Material:
        raise NotImplementedError


class PointDefectBuilder(ConvertGeneratedItemsPymatgenStructureMixin, DefectBuilder):
    """
    Builder class for generating point defects.
    """

    _BuildParametersType = PointDefectBuilderParameters
    _DefaultBuildParameters = PointDefectBuilderParameters()
    _GeneratedItemType: PymatgenStructure = PymatgenStructure
    _ConfigurationType = PointDefectConfiguration
    _generator: Callable

    def _get_species(self, configuration: BaseBuilder._ConfigurationType):
        crystal_elements = configuration.crystal.basis.elements.values
        placeholder_specie = crystal_elements[0]
        return configuration.chemical_element or placeholder_specie

    def _generate(self, configuration: BaseBuilder._ConfigurationType) -> List[_GeneratedItemType]:
        pymatgen_structure = to_pymatgen(configuration.crystal)
        pymatgen_periodic_site = PymatgenPeriodicSite(
            species=self._get_species(configuration),
            coords=configuration.coordinate,
            coords_are_cartesian=configuration.use_cartesian_coordinates,
            lattice=pymatgen_structure.lattice,
        )
        # oxi_state set to 0 to allow for 2D materials, otherwise oxi_state search takes infinite loop in pymatgen
        defect = self._generator(pymatgen_structure, pymatgen_periodic_site, oxi_state=0)
        defect_structure = defect.defect_structure.copy()
        defect_structure.remove_oxidation_states()
        return [defect_structure]

    def _update_material_name(self, material: Material, configuration: BaseBuilder._ConfigurationType) -> Material:
        updated_material = super()._update_material_name(material, configuration)
        capitalized_defect_type = configuration.defect_type.name.capitalize()
        chemical_element = configuration.chemical_element if configuration.chemical_element else ""
        new_name = f"{updated_material.name}, {capitalized_defect_type} {chemical_element} Defect"
        updated_material.name = new_name
        return updated_material


class VacancyPointDefectBuilder(PointDefectBuilder):
    _generator: PymatgenVacancy = PymatgenVacancy


class SubstitutionPointDefectBuilder(PointDefectBuilder):
    _generator: PymatgenSubstitution = PymatgenSubstitution


class InterstitialPointDefectBuilder(PointDefectBuilder):
    _generator: PymatgenInterstitial = PymatgenInterstitial


class SlabDefectBuilderParameters(BaseModel):
    auto_add_vacuum: bool = True
    vacuum_thickness: float = 5.0


class SlabDefectBuilder(DefectBuilder):
    _BuildParametersType = SlabDefectBuilderParameters
    _DefaultBuildParameters = SlabDefectBuilderParameters()

    def create_material_with_additional_layers(
        self, material: Material, added_thickness: Union[int, float] = 1
    ) -> Material:
        """
        Adds a number of layers to the material.

        Args:
            material: The original material.
            added_thickness: The thickness to add.

        Returns:
            A new Material instance with the added layers.
        """
        if isinstance(added_thickness, int):
            return self.create_material_with_additional_layers_int(material, added_thickness)
        elif isinstance(added_thickness, float):
            return self.create_material_with_additional_layers_float(material, added_thickness)
        else:
            raise TypeError("added_thickness must be an integer or float for this method.")

    def create_material_with_additional_layers_int(self, material: Material, added_thickness: int = 1) -> Material:
        """
        Adds an integer number of layers to the material.

        Args:
            material: The original material.
            added_thickness: The number of whole layers to add.

        Returns:
            A new Material instance with the added layers.
        """

        new_material = material.clone()
        termination = Termination.from_string(new_material.metadata.get("build").get("termination"))
        build_config = new_material.metadata.get("build").get("configuration")

        if build_config["type"] != "SlabConfiguration":
            raise ValueError("Material is not a slab.")
        build_config.pop("type")
        build_config["thickness"] = build_config["thickness"] + added_thickness

        new_slab_config = SlabConfiguration(**build_config)
        material_with_additional_layers = create_slab(new_slab_config, termination)

        return material_with_additional_layers

    def create_material_with_additional_layers_float(
        self, material: Material, added_thickness: float = 1.0
    ) -> Material:
        """
        Adds a fractional number of layers to the material.

        Args:
            material: The original material.
            added_thickness: The fractional thickness to add.

        Returns:
            A new Material instance with the fractional layer added.
        """
        whole_layers = int(added_thickness)
        fractional_part = added_thickness - whole_layers

        if whole_layers > 0:
            material_with_additional_layers = self.create_material_with_additional_layers_int(material, whole_layers)
        else:
            material_with_additional_layers = material.clone()

        if fractional_part > 0.0:
            material_with_additional_layers = self.add_fractional_layer(
                material_with_additional_layers, whole_layers, fractional_part
            )

        return material_with_additional_layers

    def add_fractional_layer(
        self,
        material: Material,
        whole_layers: int,
        fractional_thickness: float,
    ) -> Material:
        """
        Adds a fractional layer to the material.

        Args:
            material: The original material.
            fractional_thickness: The fractional thickness to add.

        Returns:
            A new Material instance with the fractional layer added.
        """
        material_with_additional_layers = self.create_material_with_additional_layers_int(material, 1)
        new_c = material_with_additional_layers.lattice.c
        layer_height = (new_c - material.lattice.c) / (whole_layers + 1)
        original_max_z = get_atomic_coordinates_extremum(material, "max", "z", use_cartesian_coordinates=True)
        added_layers_max_z = original_max_z + (whole_layers + fractional_thickness) * layer_height
        added_layers_max_z_crystal = material_with_additional_layers.basis.cell.convert_point_to_crystal(
            [0, 0, added_layers_max_z]
        )[2]

        material_with_additional_layers = filter_by_box(
            material=material_with_additional_layers,
            max_coordinate=[1, 1, added_layers_max_z_crystal],
        )

        return material_with_additional_layers

    def merge_slab_and_defect(self, material: Material, isolated_defect: Material) -> Material:
        new_vacuum = isolated_defect.lattice.c - material.lattice.c
        new_material = add_vacuum(material, new_vacuum)
        new_material.to_cartesian()
        new_material = merge_materials(
            materials=[isolated_defect, new_material],
            material_name=material.name,
            merge_dangerously=True,
        )
        new_material.to_crystal()
        if self.build_parameters.auto_add_vacuum and get_atomic_coordinates_extremum(new_material, "max", "z") > 1:
            new_material = add_vacuum(new_material, self.build_parameters.vacuum_thickness)
        return new_material


class AdatomSlabDefectBuilder(SlabDefectBuilder):
    _ConfigurationType: type(AdatomSlabPointDefectConfiguration) = AdatomSlabPointDefectConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    def create_adatom(
        self,
        material: Material,
        chemical_element: str = "Si",
        position_on_surface: Optional[List[float]] = None,
        distance_z: float = 2.0,
    ) -> Material:
        """
        Create an adatom at the specified position on the surface of the material.

        Args:
            material: The material to add the adatom to.
            chemical_element: The chemical element of the adatom.
            position_on_surface: The position on the surface of the material.
            distance_z: The distance of the adatom from the surface.

        Returns:
            The material with the adatom added.
        """
        if position_on_surface is None:
            position_on_surface = [0.5, 0.5]
        new_material = material.clone()
        new_basis = new_material.basis
        adatom_coordinate = self._calculate_coordinate_from_position_and_distance(
            material, position_on_surface, distance_z
        )
        new_basis.add_atom(chemical_element, adatom_coordinate)
        new_material.basis = new_basis
        return new_material

    def _calculate_coordinate_from_position_and_distance(
        self, material: Material, position_on_surface: List[float], distance_z: float
    ) -> List[float]:
        max_z_id = get_local_extremum_atom_index(
            material, position_on_surface, "max", vicinity=3.0, use_cartesian_coordinates=False
        )
        max_z = material.basis.coordinates.get_element_value_by_index(max_z_id)[2]
        distance_in_crystal_units = distance_z / material.lattice.c
        return [position_on_surface[0], position_on_surface[1], max_z + distance_in_crystal_units]

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        updated_material = super()._update_material_name(material, configuration)
        adatom_element = configuration.chemical_element
        new_name = f"{updated_material.name}, Adatom {adatom_element} Defect"
        updated_material.name = new_name
        return updated_material

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        return [
            self.create_adatom(
                material=configuration.crystal,
                chemical_element=configuration.chemical_element,
                position_on_surface=configuration.position_on_surface,
                distance_z=configuration.distance_z,
            )
        ]


class EquidistantAdatomSlabDefectBuilder(AdatomSlabDefectBuilder):
    def create_adatom(
        self,
        material: Material,
        chemical_element: str = "Si",
        position_on_surface: Optional[List[float]] = None,
        distance_z: float = 2.0,
    ) -> Material:
        """
        Create an adatom with an equidistant XY position among the nearest neighbors
        at the given distance from the surface.

        Args:
            material: The material to add the adatom to.
            chemical_element: The chemical element of the adatom.
            position_on_surface: The position on the surface of the material.
            distance_z: The distance of the adatom from the surface.

        Returns:
            The material with the adatom added.
        """
        if position_on_surface is None:
            position_on_surface = [0.5, 0.5]
        equidistant_position = self.get_equidistant_position(material, position_on_surface, distance_z)
        new_material = material.clone()
        if equidistant_position[2] > 1:
            if self.build_parameters.auto_add_vacuum:
                new_material = add_vacuum(material, self.build_parameters.vacuum_thickness)
                equidistant_position = self.get_equidistant_position(new_material, position_on_surface, distance_z)
            else:
                raise ValueError("Not enough vacuum space to place the adatom.")

        return super().create_adatom(new_material, chemical_element, equidistant_position, distance_z)

    def get_equidistant_position(
        self, material: Material, position_on_surface: List[float], distance_z: float = 2.0
    ) -> List[float]:
        adatom_coordinate = self._calculate_coordinate_from_position_and_distance(
            material, position_on_surface, distance_z
        )
        # We need to find the neighboring atoms with pbc by looking at the central unit cell in 3x3x3 supercell
        scaling_factor = [3, 3, 1]
        translation_vector = [1 / 3, 1 / 3, 0]
        supercell_material = create_supercell(material, scaling_factor=scaling_factor)
        adatom_coordinate_in_supercell = transform_coordinate_to_supercell(
            coordinate=adatom_coordinate, scaling_factor=scaling_factor, translation_vector=translation_vector
        )

        neighboring_atoms_ids_in_supercell = get_voronoi_nearest_neighbors_atom_indices(
            material=supercell_material, coordinate=adatom_coordinate_in_supercell
        )
        if neighboring_atoms_ids_in_supercell is None:
            raise ValueError("No neighboring atoms found. Try reducing the distance_z.")

        isolated_neighboring_atoms_basis = supercell_material.basis.copy()
        isolated_neighboring_atoms_basis.coordinates.filter_by_ids(neighboring_atoms_ids_in_supercell)
        equidistant_coordinate_in_supercell = get_center_of_coordinates(
            isolated_neighboring_atoms_basis.coordinates.values
        )
        equidistant_coordinate_in_supercell[2] = adatom_coordinate[2]

        return transform_coordinate_to_supercell(
            equidistant_coordinate_in_supercell, scaling_factor, translation_vector, reverse=True
        )


class CrystalSiteAdatomSlabDefectBuilder(AdatomSlabDefectBuilder):
    def calculate_approximate_adatom_coordinate(
        self, material: Material, position_on_surface: List[float], distance_z: float
    ) -> List[float]:
        approximate_adatom_coordinate = self._calculate_coordinate_from_position_and_distance(
            material, position_on_surface, distance_z
        )
        approximate_adatom_coordinate_cartesian = material.basis.cell.convert_point_to_cartesian(
            approximate_adatom_coordinate
        )
        return approximate_adatom_coordinate_cartesian

    def create_isolated_adatom(
        self,
        original_material: Material,
        configuration: AdatomSlabPointDefectConfiguration,
    ) -> Material:
        material: Material = configuration.crystal
        approximate_adatom_coordinate_cartesian: List[float] = self.calculate_approximate_adatom_coordinate(
            original_material, configuration.position_on_surface, configuration.distance_z
        )
        chemical_element: Optional[str] = configuration.chemical_element or None
        if chemical_element is None:
            closest_site_id = get_closest_site_id_from_coordinate(material, approximate_adatom_coordinate_cartesian)
        else:
            closest_site_id = get_closest_site_id_from_coordinate_and_element(
                material, approximate_adatom_coordinate_cartesian, chemical_element
            )
        only_adatom_material = filter_by_ids(material, [closest_site_id])
        return only_adatom_material

    def create_adatom(
        self,
        material: Material,
        chemical_element: Optional[str] = None,
        position_on_surface: Optional[List[float]] = None,
        distance_z: float = 0,
    ) -> Material:
        """
        Create an adatom at the crystal site closest to the specified position on the surface of the material.

        Args:
            material: The material to add the adatom to.
            chemical_element: The chemical element of the adatom.
            position_on_surface: The position on the surface of the material.
            distance_z: The distance of the adatom from the surface.

        Returns:
            The material with the adatom added.
        """
        if position_on_surface is None:
            position_on_surface = [0.5, 0.5]

        new_material = material.clone()
        material_with_additional_layer = self.create_material_with_additional_layers(new_material)
        material_with_additional_layer.to_cartesian()

        only_adatom_material = self.create_isolated_adatom(
            original_material=new_material,
            configuration=AdatomSlabPointDefectConfiguration(
                crystal=material_with_additional_layer,
                chemical_element=chemical_element,
                position_on_surface=position_on_surface,
                distance_z=distance_z,
            ),
        )

        return self.merge_slab_and_defect(new_material, only_adatom_material)


class DefectPairBuilder(DefectBuilder):
    def create_defect_pair(
        self,
        primary_defect_configuration: Union[PointDefectConfiguration, AdatomSlabPointDefectConfiguration],
        secondary_defect_configuration: Union[PointDefectConfiguration, AdatomSlabPointDefectConfiguration],
    ) -> Material:
        """
        Create a pair of point defects in the material.

        Args:
            primary_defect_configuration: The configuration of the first defect.
            secondary_defect_configuration: The configuration of the second defect.

        Returns:
            Material: The material with both defects added.
        """
        primary_material = self.create_isolated_defect(primary_defect_configuration)
        # Remove metadata to allow for independent defect creation
        if hasattr(primary_defect_configuration.crystal.metadata, "build"):
            primary_material.metadata["build"] = primary_defect_configuration.crystal.metadata["build"]
        primary_material.name = primary_defect_configuration.crystal.name
        secondary_defect_configuration.crystal = primary_material
        secondary_material = self.create_isolated_defect(secondary_defect_configuration)

        return secondary_material


class PointDefectPairBuilder(PointDefectBuilder, DefectPairBuilder):
    _ConfigurationType: type(PointDefectPairConfiguration) = PointDefectPairConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    def create_isolated_defect(self, defect_configuration: PointDefectConfiguration) -> Material:
        key = defect_configuration.defect_type.value
        if hasattr(defect_configuration, "placement_method") and defect_configuration.placement_method is not None:
            key += f":{defect_configuration.placement_method.name}".lower()
        builder_class = DefectBuilderFactory.get_class_by_name(key)
        defect_builder = builder_class()
        return defect_builder.get_material(defect_configuration)

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        return [
            self.create_defect_pair(
                primary_defect_configuration=configuration.primary_defect_configuration,
                secondary_defect_configuration=configuration.secondary_defect_configuration,
            )
        ]

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        updated_material = super()._update_material_name(material, configuration)
        name_1 = configuration.primary_defect_configuration.defect_type.name.capitalize()
        name_2 = configuration.secondary_defect_configuration.defect_type.name.capitalize()
        new_name = f"{updated_material.name}, {name_1} and {name_2} Defect Pair"
        updated_material.name = new_name
        return updated_material


class IslandSlabDefectBuilder(SlabDefectBuilder):
    _ConfigurationType: type(IslandSlabDefectConfiguration) = IslandSlabDefectConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    @staticmethod
    def _default_condition(coordinate: List[float]):
        return True

    def create_island(
        self,
        material: Material,
        condition: Optional[Callable[[List[float]], bool]] = None,
        thickness: int = 1,
        use_cartesian_coordinates: bool = False,
    ) -> Material:
        """
        Create an island at the specified position on the surface of the material.

        Args:
            material: The material to add the island to.
            condition: The condition on coordinates to describe the island.
            thickness: The thickness of the island in layers.
            use_cartesian_coordinates: Whether to use Cartesian coordinates for the condition.
        Returns:
            The material with the island added.
        """
        new_material = material.clone()
        original_max_z = get_atomic_coordinates_extremum(new_material, use_cartesian_coordinates=True)
        material_with_additional_layers = self.create_material_with_additional_layers(new_material, thickness)
        added_layers_max_z = get_atomic_coordinates_extremum(
            material_with_additional_layers, use_cartesian_coordinates=True
        )
        if condition is None:
            condition = self._default_condition

        atoms_within_island = filter_by_condition_on_coordinates(
            material=material_with_additional_layers,
            condition=condition,
            use_cartesian_coordinates=use_cartesian_coordinates,
        )
        # Filter atoms in the added layers between the original and added layers
        island_material = filter_by_box(
            material=atoms_within_island,
            min_coordinate=[0, 0, original_max_z],
            max_coordinate=[material.lattice.a, material.lattice.b, added_layers_max_z],
            use_cartesian_coordinates=True,
        )

        return self.merge_slab_and_defect(island_material, new_material)

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        condition_callable = configuration.condition.condition
        return [
            self.create_island(
                material=configuration.crystal,
                condition=condition_callable,
                thickness=configuration.number_of_added_layers,
                use_cartesian_coordinates=configuration.use_cartesian_coordinates,
            )
        ]


class TerraceSlabDefectBuilder(SlabDefectBuilder):
    _ConfigurationType: type(TerraceSlabDefectConfiguration) = TerraceSlabDefectConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    def _calculate_cut_direction_vector(self, material: Material, cut_direction: List[int]):
        """
        Calculate the normalized cut direction vector in Cartesian coordinates.

        Args:
            material: The material to get the lattice vectors from.
            cut_direction: The direction of the cut in lattice directions.

        Returns:
            The normalized cut direction vector in Cartesian coordinates.
        """
        np_cut_direction = np.array(cut_direction)
        direction_vector = np.dot(np.array(material.basis.cell.vectors_as_array), np_cut_direction)
        normalized_direction_vector = direction_vector / np.linalg.norm(direction_vector)
        return normalized_direction_vector

    def _calculate_height_cartesian(self, material: Material, new_material: Material):
        """
        Calculate the height of the added layers in Cartesian coordinates.

        Args:
            material: The original material.
            new_material: The material with the added layers.

        Returns:
            The height of the added layers in Cartesian coordinates.
        """
        original_max_z = get_atomic_coordinates_extremum(material, use_cartesian_coordinates=True)
        added_layers_max_z = get_atomic_coordinates_extremum(new_material, use_cartesian_coordinates=True)
        height_cartesian = added_layers_max_z - original_max_z
        return height_cartesian

    def _calculate_rotation_parameters(
        self, original_material: Material, new_material: Material, normalized_direction_vector: List[float]
    ):
        """
        Calculate the necessary rotation angle and axis.

        Args:
            original_material: The original material.
            new_material: The material with the added layers.
            normalized_direction_vector: The normalized cut direction vector in Cartesian coordinates.

        Returns:
            The rotation angle, normalized rotation axis, and delta length.
        """
        height_cartesian = self._calculate_height_cartesian(original_material, new_material)
        cut_direction_xy_proj_cart = np.linalg.norm(
            np.dot(np.array(new_material.basis.cell.vectors_as_array), normalized_direction_vector)
        )
        # Slope of the terrace along the cut direction
        hypotenuse = np.linalg.norm([height_cartesian, cut_direction_xy_proj_cart])
        angle = np.arctan(height_cartesian / cut_direction_xy_proj_cart) * 180 / np.pi
        normalized_rotation_axis = np.cross(normalized_direction_vector, [0, 0, 1]).tolist()
        delta_length = hypotenuse - cut_direction_xy_proj_cart
        return angle, normalized_rotation_axis, delta_length

    def _increase_lattice_size(
        self, material: Material, length_increase: float, direction_of_increase: List[float]
    ) -> Material:
        """
        Increase the lattice size in a specific direction.

        When the material is rotated to maintain periodic boundary conditions (PBC),
        the periodicity in the X and Y directions changes.
        Therefore, the lattice size must be increased to fit the new structure dimensions.

        If the terrace plane is normal to the Z direction, it becomes larger than the previous XY plane of the material
        because it forms a hypotenuse between PBC points.
        This method adjusts the lattice vectors to accommodate this change.

        Args:
            material: The material to increase the lattice size of.
            length_increase: The increase in length.
            direction_of_increase: The direction of the increase.

        Returns:
            The material with the increased lattice size.
        """
        vector_a, vector_b = np.array(material.basis.cell.vector1), np.array(material.basis.cell.vector2)
        norm_a, norm_b = np.linalg.norm(vector_a), np.linalg.norm(vector_b)

        delta_a_cart = np.dot(vector_a, np.array(direction_of_increase)) * length_increase / norm_a
        delta_b_cart = np.dot(vector_b, np.array(direction_of_increase)) * length_increase / norm_b

        scaling_matrix = np.eye(3)
        scaling_matrix[0, 0] += delta_a_cart / norm_a
        scaling_matrix[1, 1] += delta_b_cart / norm_b

        cart_basis = material.basis.copy()
        cart_basis.to_cartesian()
        cart_basis.cell.scale_by_matrix(scaling_matrix)
        material.basis = cart_basis

        new_lattice = material.lattice.clone()
        new_lattice.a = np.linalg.norm(cart_basis.cell.vector1)
        new_lattice.b = np.linalg.norm(cart_basis.cell.vector2)
        material.lattice = new_lattice
        return material

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        new_material = super()._update_material_name(material, configuration)
        new_name = (
            f"{new_material.name}, {configuration.number_of_added_layers}-step Terrace {configuration.cut_direction}"
        )
        new_material.name = new_name
        return new_material

    def create_terrace(
        self,
        material: Material,
        cut_direction: Optional[List[int]] = None,
        pivot_coordinate: Optional[List[float]] = None,
        number_of_added_layers: int = 1,
        use_cartesian_coordinates: bool = False,
        rotate_to_match_pbc: bool = True,
    ) -> Material:
        """
        Create a terrace at the specified position on the surface of the material.

        Args:
            material: The material to add the terrace to.
            cut_direction: The direction of the cut in lattice directions.
            pivot_coordinate: The center position of the terrace.
            number_of_added_layers: The number of added layers to the slab which will form the terrace
            use_cartesian_coordinates: Whether to use Cartesian coordinates for the center position.
            rotate_to_match_pbc: Whether to rotate the material to match the periodic boundary conditions.
        Returns:
            The material with the terrace added.
        """
        if cut_direction is None:
            cut_direction = [0, 0, 1]
        if pivot_coordinate is None:
            pivot_coordinate = [0.5, 0.5, 0.5]

        new_material = material.clone()
        material_with_additional_layers = self.create_material_with_additional_layers(
            new_material, number_of_added_layers
        )

        normalized_direction_vector = self._calculate_cut_direction_vector(material, cut_direction)
        condition = CoordinateCondition.PlaneCoordinateCondition(
            plane_normal=normalized_direction_vector,
            plane_point_coordinate=pivot_coordinate,
        ).condition
        atoms_within_terrace = filter_by_condition_on_coordinates(
            material=material_with_additional_layers,
            condition=condition,
            use_cartesian_coordinates=use_cartesian_coordinates,
        )
        merged_material = self.merge_slab_and_defect(new_material, atoms_within_terrace)

        angle, normalized_rotation_axis, delta_length = self._calculate_rotation_parameters(
            material, merged_material, normalized_direction_vector
        )
        result_material = translate_to_z_level(merged_material, "center")

        if rotate_to_match_pbc:
            adjusted_material = self._increase_lattice_size(result_material, delta_length, normalized_direction_vector)
            result_material = rotate(material=adjusted_material, axis=normalized_rotation_axis, angle=angle)
        return result_material

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        return [
            self.create_terrace(
                material=configuration.crystal,
                cut_direction=configuration.cut_direction,
                pivot_coordinate=configuration.pivot_coordinate,
                number_of_added_layers=configuration.number_of_added_layers,
                use_cartesian_coordinates=configuration.use_cartesian_coordinates,
            )
        ]
