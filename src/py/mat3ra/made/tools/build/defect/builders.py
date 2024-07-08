from typing import List, Callable, Optional

from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.modify import add_vacuum
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
from ...analyze import get_nearest_neighbors_atom_indices, get_atomic_coordinates_extremum
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
        adatom_position = self._calculate_position_from_2d(material, position_on_surface, distance_z)
        new_basis.add_atom(chemical_element, adatom_position)
        new_material.basis = new_basis
        return [new_material]

    def _calculate_position_from_2d(
        self, material: Material, position_on_surface: List[float], distance_z: float
    ) -> List[float]:
        max_z = get_atomic_coordinates_extremum(material)
        distance_z = distance_z
        distance_in_crystal_units = distance_z / material.lattice.c
        position = position_on_surface.copy()
        position = position[:2]
        position.append(max_z + distance_in_crystal_units)
        return position

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
        new_basis = material.basis
        adatom_position = self._calculate_position_from_2d(material, position_on_surface, distance_z)
        neighboring_atoms_ids = get_nearest_neighbors_atom_indices(material, adatom_position)
        # We need to check if neighboring atoms number is the same in pbc
        supercell_material = create_supercell(material, [[3, 0, 0], [0, 3, 0], [0, 0, 1]])
        # Move the coordinate to the central unit cell of the supercell (crystal coordinates)
        supercell_adatom_position = [1 / 3 + adatom_position[0] / 3, 1 / 3 + adatom_position[1] / 3, adatom_position[2]]
        supercell_neighboring_atoms_ids = get_nearest_neighbors_atom_indices(
            supercell_material, supercell_adatom_position
        )
        if neighboring_atoms_ids is None or supercell_neighboring_atoms_ids is None:
            raise ValueError("No neighboring atoms found. Try reducing the distance_z.")
        if len(supercell_neighboring_atoms_ids) != len(neighboring_atoms_ids):
            raise ValueError("Number of neighboring atoms is not the same in PBC. Try increasing the supercell size.")
        neighboring_atoms_coordinates = [new_basis.coordinates.values[atom_id] for atom_id in neighboring_atoms_ids]

        equidistant_position = get_center_of_coordinates(neighboring_atoms_coordinates)
        equidistant_position[2] = adatom_position[2]

        return equidistant_position
