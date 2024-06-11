from typing import List, Callable

from pydantic import BaseModel
from pymatgen.analysis.defects.core import (
    Substitution as PymatgenSubstitution,
    Vacancy as PymatgenVacancy,
    Interstitial as PymatgenInterstitial,
)
from pymatgen.core import PeriodicSite as PymatgenPeriodicSite

from ...build import BaseBuilder
from ...convert import PymatgenStructure, to_pymatgen
from ..mixins import ConvertGeneratedItemsPymatgenStructureMixin
from .configuration import PointDefectConfiguration
from mat3ra.made.utils import get_array_with_id_value_element_value_by_index


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

    def _get_species(self, configuration: PointDefectConfiguration):
        crystal_elements = configuration.crystal.basis["elements"]
        placeholder_specie = get_array_with_id_value_element_value_by_index(crystal_elements, 0)
        return configuration.chemical_element or placeholder_specie

    def _generate(self, configuration: PointDefectConfiguration) -> List[_GeneratedItemType]:
        pymatgen_structure = to_pymatgen(configuration.crystal)
        pymatgen_periodic_site = PymatgenPeriodicSite(
            species=self._get_species(configuration),
            coords=configuration.position,
            lattice=pymatgen_structure.lattice,
        )
        defect = self._generator(pymatgen_structure, pymatgen_periodic_site)
        return [defect.defect_structure]


class VacancyBuilder(PointDefectBuilder):
    _generator: PymatgenVacancy = PymatgenVacancy


class SubstitutionBuilder(PointDefectBuilder):
    _generator: PymatgenSubstitution = PymatgenSubstitution


class InterstitialBuilder(PointDefectBuilder):
    _generator: PymatgenInterstitial = PymatgenInterstitial