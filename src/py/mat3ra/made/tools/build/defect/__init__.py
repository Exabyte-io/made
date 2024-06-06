from enum import Enum
from typing import Optional, List, Union, Any

from pydantic import BaseModel

from ...build import BaseBuilder
from ...convert import PymatgenStructure, to_pymatgen, from_pymatgen
from pymatgen.analysis.defects.core import Substitution, Vacancy, Interstitial
from pymatgen.core import PeriodicSite


class PointDefectTypeEnum(str, Enum):
    VACANCY = "vacancy"
    SUBSTITUTION = "substitution"
    INTERSTITIAL = "interstitial"


class BaseDefectConfiguration(BaseModel):
    # TODO: fix arbitrary_types_allowed error and set Material class type
    material: Any
    # defect_type type can be an Enum for a specific defect class (for point defect, 2d defect, etc.)
    defect_type: Union[PointDefectTypeEnum, None] = None


class PointDefectConfiguration(BaseDefectConfiguration):
    defect_type: PointDefectTypeEnum = PointDefectTypeEnum.VACANCY
    # TODO: should come from Enum of elements
    specie: Optional[str] = None


class PointDefectBuilderParameters(BaseModel):
    target_site: int = 0
    # TODO: import coordinate type from ESSE
    position: Optional[List[int]] = None
    center_defect: bool = False


class PointDefectBuilder(BaseBuilder):
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
        pymatgen_periodic_site = PeriodicSite(
            species=pymatgen_structure.species[site_index],
            coords=pymatgen_site.frac_coords,
            lattice=pymatgen_structure.lattice,
        )

        if configuration.defect_type == PointDefectTypeEnum.VACANCY:
            defect = Vacancy(pymatgen_structure, pymatgen_periodic_site)
        elif configuration.defect_type == PointDefectTypeEnum.SUBSTITUTION:
            defect = Substitution(pymatgen_structure, pymatgen_periodic_site, configuration.specie)
        elif configuration.defect_type == PointDefectTypeEnum.INTERSTITIAL:
            defect = Interstitial(pymatgen_structure, pymatgen_periodic_site)
        else:
            raise ValueError(f"Unknown defect type: {configuration.defect_type}")

        return [defect.centered_defect_structure if self.build_parameters.center_defect else defect.defect_structure]

    @staticmethod
    def _convert_generated_item(item: _GeneratedItemType):
        return from_pymatgen(item)
