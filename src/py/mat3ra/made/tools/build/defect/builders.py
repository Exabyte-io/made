from typing import List

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
from .configuration import PointDefectConfiguration


class PointDefectBuilderParameters(BaseModel):
    target_site: int = 0
    center_defect: bool = False


class PointDefectBuilder(ConvertGeneratedItemsPymatgenStructureMixin, BaseBuilder):
    """
    Builder class for generating point defects.
    """

    _BuildParametersType = PointDefectBuilderParameters
    _DefaultBuildParameters = PointDefectBuilderParameters()
    _GeneratedItemType: PymatgenStructure = PymatgenStructure
    _ConfigurationType = PointDefectConfiguration

    def _generate(self, configuration: PointDefectConfiguration) -> List[_GeneratedItemType]:
        site_index = self.build_parameters.target_site
        pymatgen_structure = to_pymatgen(configuration.material)
        pymatgen_site = pymatgen_structure[site_index]
        interstitial_frac_coords_shift = configuration.position_shift if configuration.position_shift else [0, 0, 0]
        pymatgen_periodic_site = PymatgenPeriodicSite(
            species=configuration.specie if configuration.specie else pymatgen_site.specie,
            coords=pymatgen_site.frac_coords + interstitial_frac_coords_shift,
            lattice=pymatgen_structure.lattice,
        )

        if configuration.defect_type == PointDefectTypeEnum.VACANCY:
            defect = PymatgenVacancy(structure=pymatgen_structure, site=pymatgen_periodic_site)
        elif configuration.defect_type == PointDefectTypeEnum.SUBSTITUTION:
            defect = PymatgenSubstitution(structure=pymatgen_structure, site=pymatgen_periodic_site)
        elif configuration.defect_type == PointDefectTypeEnum.INTERSTITIAL:
            defect = PymatgenInterstitial(structure=pymatgen_structure, site=pymatgen_periodic_site)
        else:
            raise ValueError(f"Unknown defect type: {configuration.defect_type}")

        return [defect.centered_defect_structure if self.build_parameters.center_defect else defect.defect_structure]
