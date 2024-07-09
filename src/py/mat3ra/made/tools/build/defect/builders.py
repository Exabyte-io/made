from typing import List, Callable, Optional

from mat3ra.made.tools.build.slab import SlabConfiguration, create_slab, Termination
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.build.utils import merge_materials
from mat3ra.made.tools.modify import add_vacuum, filter_material_by_ids
from pydantic import BaseModel
from mat3ra.made.material import Material

from ...third_party import (
    PymatgenStructure,
    PymatgenPeriodicSite,
    PymatgenVacancy,
    PymatgenSubstitution,
    PymatgenInterstitial,
)
from ...build import BaseBuilder
from ...convert import to_pymatgen
from ...analyze import (
    get_nearest_neighbors_atom_indices,
    get_atomic_coordinates_extremum,
    get_closest_site_id_from_position,
    get_closest_site_id_from_position_and_element,
)
from ....utils import get_center_of_coordinates
from ..mixins import ConvertGeneratedItemsPymatgenStructureMixin
from .configuration import PointDefectConfiguration, AdatomSlabDefectConfiguration


class PointDefectBuilderParameters(BaseModel):
    center_defect: bool = False


class PointDefectBuilder(ConvertGeneratedItemsPymatgenStructureMixin, BaseBuilder):
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
            coords=configuration.position,
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


class SlabDefectBuilder(BaseBuilder):
    _BuildParametersType = SlabDefectBuilderParameters
    _DefaultBuildParameters = SlabDefectBuilderParameters()


class AdatomSlabDefectBuilder(SlabDefectBuilder):
    _ConfigurationType: type(AdatomSlabDefectConfiguration) = AdatomSlabDefectConfiguration  # type: ignore
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
        coordinate = position_on_surface.copy()
        coordinate = coordinate[:2]
        coordinate.append(max_z + distance_in_crystal_units)
        return coordinate

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
        new_basis = material.basis.copy()
        adatom_coordinate = self._calculate_coordinate_from_position_and_distance(
            material, position_on_surface, distance_z
        )
        neighboring_atoms_ids = get_nearest_neighbors_atom_indices(material, adatom_coordinate)
        # We need to check if neighboring atoms number is the same in pbc
        supercell_material = create_supercell(material, [[3, 0, 0], [0, 3, 0], [0, 0, 1]])
        # Move the coordinate to the central unit cell of the supercell (crystal coordinates)
        supercell_adatom_coordinate = [
            1 / 3 + adatom_coordinate[0] / 3,
            1 / 3 + adatom_coordinate[1] / 3,
            adatom_coordinate[2],
        ]
        supercell_neighboring_atoms_ids = get_nearest_neighbors_atom_indices(
            supercell_material, supercell_adatom_coordinate
        )
        if neighboring_atoms_ids is None or supercell_neighboring_atoms_ids is None:
            raise ValueError("No neighboring atoms found. Try reducing the distance_z.")
        if len(supercell_neighboring_atoms_ids) != len(neighboring_atoms_ids):
            raise ValueError("Number of neighboring atoms is not the same in PBC. Try increasing the supercell size.")
        print("lattice", material.lattice)
        print("coordinates pre ", new_basis.coordinates.to_array_of_values_with_ids())
        new_basis.coordinates.filter_by_ids(neighboring_atoms_ids)
        neighboring_atoms_coordinates = new_basis.coordinates.values
        print("coordinates post", new_basis.coordinates.to_array_of_values_with_ids())
        equidistant_coordinate = get_center_of_coordinates(neighboring_atoms_coordinates)
        equidistant_coordinate[2] = adatom_coordinate[2]

        return equidistant_coordinate


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

    def create_material_with_additional_layer(self, material: Material) -> Material:
        termination = Termination.from_string(material.metadata.get("build").get("termination"))
        build_config = material.metadata.get("build").get("configuration")
        if build_config["type"] != "SlabConfiguration":
            raise ValueError("Material is not a slab.")
        build_config.pop("type")
        build_config["thickness"] = build_config["thickness"] + 1
        new_slab_config = SlabConfiguration(**build_config)
        material_with_additional_layer = create_slab(new_slab_config, termination)

        cartesian_basis = material_with_additional_layer.basis
        cartesian_basis.to_cartesian()
        material_with_additional_layer.basis = cartesian_basis

        return material_with_additional_layer

    def create_only_adatom_material(
        self,
        material_with_additional_layer: Material,
        approximate_adatom_coordinate_cartesian: List[float],
        chemical_element: Optional[str] = None,
    ) -> Material:
        if chemical_element is None:
            closest_site_id = get_closest_site_id_from_position(
                material_with_additional_layer, approximate_adatom_coordinate_cartesian
            )
        else:
            closest_site_id = get_closest_site_id_from_position_and_element(
                material_with_additional_layer, approximate_adatom_coordinate_cartesian, chemical_element
            )
        only_adatom_material = filter_material_by_ids(material_with_additional_layer, [closest_site_id])
        return only_adatom_material

    def merge_material_with_adatom(self, material: Material, only_adatom_material: Material) -> Material:
        new_vacuum = only_adatom_material.lattice.c - material.lattice.c
        new_material = add_vacuum(material, new_vacuum)
        cartesian_basis = new_material.basis
        cartesian_basis.to_cartesian()
        new_material.basis = cartesian_basis
        new_material = merge_materials(
            [only_adatom_material, new_material],
            merge_dangerously=True,
        )
        return new_material

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

        material_with_additional_layer = self.create_material_with_additional_layer(new_material)

        only_adatom_material = self.create_only_adatom_material(
            material_with_additional_layer, approximate_adatom_coordinate_cartesian, chemical_element
        )

        return [self.merge_material_with_adatom(new_material, only_adatom_material)]
