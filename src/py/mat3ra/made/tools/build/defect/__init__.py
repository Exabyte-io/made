from enum import Enum
from typing import Optional, List, Literal, Union

from pydantic import BaseModel

from src.py.mat3ra.made.material import Material
from ...build import BaseBuilder
from ...convert import PymatgenStructure, to_pymatgen
from pymatgen.analysis.defects.core import Substitution, Vacancy, Interstitial
from pymatgen.core import PeriodicSite, Species, Structure


class PointDefectTypeEnum(str, Enum):
    VACANCY = "vacancy"
    SUBSTITUTION = "substitution"
    INTERSTITIAL = "interstitial"


class BaseDefectConfiguration(BaseModel):
    material: Material
    # defect_type type can be an Enum for a specific defect class (for point defect, 2d defect, etc.)
    defect_type: Union[PointDefectTypeEnum, None] = None


class PointDefectConfiguration(BaseDefectConfiguration):
    defect_type: PointDefectTypeEnum = PointDefectTypeEnum.VACANCY
    # TODO: should come from Enum of elements
    specie: Optional[str] = None
    min_distance: Optional[float] = 0.0


class PointDefectBuilderParameters(BaseModel):
    target_site: int = 0
    # TODO: import coordinate type from ESSE
    position: Optional[List[int]] = None
    center_defect: bool = False
    resize_cell_matrix: Union[List[int], Literal["auto"]] = "auto"


class PointDefectBuilder(BaseBuilder):
    """
    Builder class for generating point defects.
    """

    _BuildParametersType = PointDefectBuilderParameters
    _DefaultBuildParameters = PointDefectBuilderParameters()
    _GeneratedItemType = PymatgenStructure
    _ConfigurationType = PointDefectConfiguration
    _SelectorParametersType = None
    _PostProcessParametersType = None

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

        return [defect.defect_structure]
