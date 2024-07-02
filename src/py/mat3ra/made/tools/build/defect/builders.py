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
from ...analyze import get_neighboring_atoms_indices


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
    add_vacuum: bool = True
    min_vacuum_thickness: float = 5.0


class SlabDefectBuilder(BaseBuilder):
    _BuildParametersType = SlabDefectBuilderParameters
    _DefaultBuildParameters = SlabDefectBuilderParameters()


class AdatomSlabDefectBuilder(SlabDefectBuilder):
    _GeneratedItemType: Material = Material

    def add_adatom(
        material: Material, chemical_element: str = "Si", position: List[float] = [0.5, 0.5, 0.5]
    ) -> Material:
        material_copy = material.clone()
        basis = material_copy.basis
        basis.add_atom(chemical_element, position)
        material_copy.basis = basis
        return [material_copy]

    _generator = add_adatom


class EquidistantAdatomSlabDefectBuilder(SlabDefectBuilder):
    _GeneratedItemType: Material = Material

    def add_adatom_equdistant(
        material: Material,
        chemical_element: str = "Si",
        approximate_position: List[float] = [0.5, 0.5, 0.5],
        distance: float = 2.0,
    ) -> Material:
        """
        Adds an atom to the material at a position that is equidistant to the nearest atoms
        (that are found within proximity to the approx position) with the specified distance.
        """
        material_copy: Material = material.clone()
        basis = material_copy.basis

        distance = distance / material_copy.lattice.c
        max_z = max([coordinate[2] for coordinate in basis.coordinates.values])

        # Move the appx position to the top of the slab
        adatom_position = approximate_position.copy()
        adatom_position[2] = max_z + distance

        neighboring_atoms_ids = get_neighboring_atoms_indices(material, adatom_position)

        # Find equidistant position from heighboring atoms with the set z
        if neighboring_atoms_ids is not None:
            neighboring_atoms_coordinates = [basis.coordinates.values[atom_id] for atom_id in neighboring_atoms_ids]

        equidistant_position = [
            sum([coordinate[i] for coordinate in neighboring_atoms_coordinates]) / len(neighboring_atoms_coordinates)
            for i in range(3)
        ]
        equidistant_position[2] = adatom_position[2]

        # Check for valid position (inside the cell)
        if equidistant_position[2] > basis.cell.vectors_as_nested_array[2][2]:
            raise ValueError("The adatom position is outside the cell.")

        basis.add_atom(chemical_element, equidistant_position)
        material_copy.basis = basis
        return [material_copy]

    _generator = add_adatom_equdistant
