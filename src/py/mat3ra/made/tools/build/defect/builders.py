from typing import List, Optional

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
from .enums import PointDefectTypeEnum
from .configuration import (
    BasePointDefectConfiguration,
    VacancyConfiguration,
    SubstitutionConfiguration,
    InterstitialConfiguration,
)


class PointDefectBuilderParameters(BaseModel):
    center_defect: bool = False


class PointDefectBuilder(ConvertGeneratedItemsPymatgenStructureMixin, BaseBuilder):
    """
    Builder class for generating point defects.
    """

    _BuildParametersType = PointDefectBuilderParameters
    _DefaultBuildParameters = PointDefectBuilderParameters()
    _GeneratedItemType: PymatgenStructure = PymatgenStructure
    pymatgen_periodic_site: PymatgenPeriodicSite = None
    #     return [defect.centered_defect_structure if self.build_parameters.center_defect else defect.defect_structure]


class VacancyBuilder(PointDefectBuilder):
    _ConfigurationType: type(VacancyConfiguration) = VacancyConfiguration  # type: ignore

    def _generate(self, configuration: _ConfigurationType) -> List[PointDefectBuilder._GeneratedItemType]:
        pymatgen_structure = to_pymatgen(configuration.crystal)
        pymatgen_site = pymatgen_structure[configuration.site_id]
        self.pymatgen_periodic_site = PymatgenPeriodicSite(
            species=pymatgen_site.specie,
            coords=pymatgen_site.frac_coords,
            lattice=pymatgen_structure.lattice,
        )
        return [PymatgenVacancy(structure=pymatgen_structure, site=self.pymatgen_periodic_site).defect_structure]


class SubstitutionBuilder(PointDefectBuilder):
    _ConfigurationType: type(SubstitutionConfiguration) = SubstitutionConfiguration  # type: ignore

    def _generate(self, configuration: _ConfigurationType) -> List[PointDefectBuilder._GeneratedItemType]:
        pymatgen_structure = to_pymatgen(configuration.crystal)
        pymatgen_site = pymatgen_structure[configuration.site_id]
        self.pymatgen_periodic_site = PymatgenPeriodicSite(
            species=configuration.element,
            coords=pymatgen_site.frac_coords,
            lattice=pymatgen_structure.lattice,
        )
        return [PymatgenSubstitution(structure=pymatgen_structure, site=self.pymatgen_periodic_site).defect_structure]


class InterstitialBuilder(PointDefectBuilder):
    _ConfigurationType: type(InterstitialConfiguration) = InterstitialConfiguration  # type: ignore

    def _generate(self, configuration: _ConfigurationType) -> List[PointDefectBuilder._GeneratedItemType]:
        pymatgen_structure = to_pymatgen(configuration.crystal)
        self.pymatgen_periodic_site = PymatgenPeriodicSite(
            species=configuration.element,
            coords=configuration.position,
            lattice=pymatgen_structure.lattice,
        )
        return [PymatgenInterstitial(structure=pymatgen_structure, site=self.pymatgen_periodic_site).defect_structure]
