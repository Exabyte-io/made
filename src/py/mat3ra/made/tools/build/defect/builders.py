from typing import List, Callable

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
from ..mixins import ConvertGeneratedItemsPymatgenStructureMixin
from .configuration import PointDefectConfiguration
from ...analyze import get_nearest_neighbors_atom_indices, get_center_of_coordinates


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
    _GeneratedItemType: Material = Material

    def add_adatom(
        self,
        material: Material,
        chemical_element: str = "Si",
        position_on_surface: List[float] = [0.5, 0.5],
        distance_z: float = 2.0,
    ) -> List[Material]:
        material_copy = material.clone()
        basis = material_copy.basis
        distance_in_crystal_units = distance_z / material_copy.lattice.c
        max_z = max([coordinate[2] for coordinate in basis.coordinates.values])
        position = position_on_surface.copy()
        position[2] = max_z + distance_in_crystal_units
        basis.add_atom(chemical_element, position)
        material_copy.basis = basis
        return [material_copy]

    _generator = add_adatom


class EquidistantAdatomSlabDefectBuilder(SlabDefectBuilder):
    _GeneratedItemType: Material = Material

    def add_adatom_equdistant(
        self,
        material: Material,
        chemical_element: str = "Si",
        approximate_position_on_surface: List[float] = [0.5, 0.5],
        distance_z: float = 2.0,
    ) -> List[Material]:
        """
        Add an atom to the material at a position that is equidistant to the nearest atoms
        (that are found within proximity to the approx position) with the specified distance.

        Args:
            material (Material): The material to add the adatom to.
            chemical_element (str): The chemical element of the adatom.
            approximate_position_on_surface (List[float]): The approximate position of the adatom on the surface.
            distance_z (float): The distance from the nearest atoms to the adatom.

        Returns:
            List[Material]: A list containing the material with the adatom added.
        """
        material_copy: Material = material.clone()
        basis = material_copy.basis
        distance_in_crystal_units = distance_z / material_copy.lattice.c
        max_z = max([coordinate[2] for coordinate in basis.coordinates.values])
        adatom_position = approximate_position_on_surface.copy()
        adatom_position[2] = max_z + distance_in_crystal_units

        neighboring_atoms_ids = get_nearest_neighbors_atom_indices(material, adatom_position)
        if not neighboring_atoms_ids:
            raise ValueError("No neighboring atoms found.")
        neighboring_atoms_coordinates = [basis.coordinates.values[atom_id] for atom_id in neighboring_atoms_ids]

        equidistant_position = get_center_of_coordinates(neighboring_atoms_coordinates)
        equidistant_position[2] = adatom_position[2]
        if equidistant_position[2] > basis.cell.vectors_as_nested_array[2][2]:
            raise ValueError("The adatom position is outside the cell.")

        basis.add_atom(chemical_element, equidistant_position)
        material_copy.basis = basis
        return [material_copy]

    _generator = add_adatom_equdistant
