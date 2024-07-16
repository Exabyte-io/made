from typing import List, Callable, Optional

from pydantic import BaseModel
from mat3ra.made.material import Material

from ...third_party import (
    PymatgenStructure,
    PymatgenPeriodicSite,
    PymatgenVacancy,
    PymatgenSubstitution,
    PymatgenInterstitial,
)

from ...modify import (
    add_vacuum,
    filter_material_by_ids,
    filter_by_box,
    filter_by_condition_on_coordinates,
)
from ...build import BaseBuilder
from ...convert import to_pymatgen
from ...analyze import (
    get_nearest_neighbors_atom_indices,
    get_atomic_coordinates_extremum,
    get_closest_site_id_from_coordinate,
    get_closest_site_id_from_coordinate_and_element,
)
from ....utils import get_center_of_coordinates
from ...utils import transform_coordinate_to_supercell
from ..utils import merge_materials
from ..slab import SlabConfiguration, create_slab, Termination
from ..supercell import create_supercell
from ..mixins import ConvertGeneratedItemsPymatgenStructureMixin
from .configuration import PointDefectConfiguration, AdatomSlabPointDefectConfiguration, IslandSlabDefectConfiguration


class PointDefectBuilderParameters(BaseModel):
    center_defect: bool = False


class DefectBuilder(BaseBuilder):
    def create_isolated_defect(self, material: Material, defect_configuration: PointDefectConfiguration) -> Material:
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
            lattice=pymatgen_structure.lattice,
        )
        defect = self._generator(pymatgen_structure, pymatgen_periodic_site)
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

    def create_material_with_additional_layers(self, material: Material, added_thickness: int = 1) -> Material:
        new_material = material.clone()
        termination = Termination.from_string(new_material.metadata.get("build").get("termination"))
        build_config = new_material.metadata.get("build").get("configuration")
        if build_config["type"] != "SlabConfiguration":
            raise ValueError("Material is not a slab.")
        build_config.pop("type")
        build_config["thickness"] = build_config["thickness"] + added_thickness
        new_slab_config = SlabConfiguration(**build_config)
        material_with_additional_layer = create_slab(new_slab_config, termination)

        return material_with_additional_layer

    def merge_slab_and_defect(self, material: Material, isolated_defect: Material) -> Material:
        new_vacuum = isolated_defect.lattice.c - material.lattice.c
        new_material = add_vacuum(material, new_vacuum)
        new_material.to_cartesian()
        new_material = merge_materials(
            materials=[isolated_defect, new_material],
            material_name=material.name,
            merge_dangerously=True,
        )
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
    ) -> List[Material]:
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
        return [new_material]

    def _calculate_coordinate_from_position_and_distance(
        self, material: Material, position_on_surface: List[float], distance_z: float
    ) -> List[float]:
        max_z = get_atomic_coordinates_extremum(material)
        distance_z = distance_z
        distance_in_crystal_units = distance_z / material.lattice.c
        return [position_on_surface[0], position_on_surface[1], max_z + distance_in_crystal_units]

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        updated_material = super()._update_material_name(material, configuration)
        adatom_element = configuration.chemical_element
        new_name = f"{updated_material.name}, Adatom {adatom_element} Defect"
        updated_material.name = new_name
        return updated_material

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        return self.create_adatom(
            material=configuration.crystal,
            chemical_element=configuration.chemical_element,
            position_on_surface=configuration.position_on_surface,
            distance_z=configuration.distance_z,
        )


class EquidistantAdatomSlabDefectBuilder(AdatomSlabDefectBuilder):
    def create_adatom(
        self,
        material: Material,
        chemical_element: str = "Si",
        position_on_surface: Optional[List[float]] = None,
        distance_z: float = 2.0,
    ) -> List[Material]:
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

        neighboring_atoms_ids_in_supercell = get_nearest_neighbors_atom_indices(
            supercell_material, adatom_coordinate_in_supercell
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

    def create_isolated_defect(
        self,
        material: Material,
        approximate_adatom_coordinate_cartesian: List[float],
        chemical_element: Optional[str] = None,
    ) -> Material:
        if chemical_element is None:
            closest_site_id = get_closest_site_id_from_coordinate(material, approximate_adatom_coordinate_cartesian)
        else:
            closest_site_id = get_closest_site_id_from_coordinate_and_element(
                material, approximate_adatom_coordinate_cartesian, chemical_element
            )
        only_adatom_material = filter_material_by_ids(material, [closest_site_id])
        return only_adatom_material

    def create_adatom(
        self,
        material: Material,
        chemical_element: Optional[str] = None,
        position_on_surface: Optional[List[float]] = None,
        distance_z: float = 0,
    ) -> List[Material]:
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

        approximate_adatom_coordinate_cartesian = self.calculate_approximate_adatom_coordinate(
            new_material, position_on_surface, distance_z
        )

        material_with_additional_layer = self.create_material_with_additional_layers(new_material)
        material_with_additional_layer.to_cartesian()

        only_adatom_material = self.create_isolated_defect(
            material_with_additional_layer, approximate_adatom_coordinate_cartesian, chemical_element
        )

        return [self.merge_slab_and_defect(new_material, only_adatom_material)]


class IslandSlabDefectBuilder(SlabDefectBuilder):
    _ConfigurationType: type(IslandSlabDefectConfiguration) = IslandSlabDefectConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    def create_island(
        self,
        material: Material,
        condition: Optional[Callable[[List[float]], bool]] = None,
        thickness: int = 1,
        use_cartesian_coordinates: bool = False,
    ) -> List[Material]:
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
        original_max_z = get_atomic_coordinates_extremum(new_material, use_cartesian_coordinates=False)
        material_with_additional_layers = self.create_material_with_additional_layers(new_material, thickness)
        added_layers_max_z = get_atomic_coordinates_extremum(material_with_additional_layers)

        if condition is None:

            def condition(coordinate: List[float]):
                return True

        atoms_within_island = filter_by_condition_on_coordinates(
            material=material_with_additional_layers,
            condition=condition,
            use_cartesian_coordinates=use_cartesian_coordinates,
        )

        # Filter atoms in the added layers
        island_material = filter_by_box(
            material=atoms_within_island,
            min_coordinate=[0, 0, original_max_z],
            max_coordinate=[1, 1, added_layers_max_z],
        )

        return [self.merge_slab_and_defect(island_material, new_material)]

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        condition_callable, _ = configuration.condition
        return self.create_island(
            material=configuration.crystal,
            condition=condition_callable,
            thickness=configuration.thickness,
            use_cartesian_coordinates=configuration.use_cartesian_coordinates,
        )
