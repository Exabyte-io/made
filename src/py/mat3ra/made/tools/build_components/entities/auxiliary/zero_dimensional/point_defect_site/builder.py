from ......build_components import MaterialWithBuildMetadata
from ..... import BaseSingleBuilder
from .configuration import PointDefectSiteConfiguration


class PointDefectSiteBuilder(BaseSingleBuilder):
    """
    Builder class for creating a material from a PointDefectSite configuration.
    """

    _ConfigurationType = PointDefectSiteConfiguration

    def _generate(self, configuration: PointDefectSiteConfiguration) -> MaterialWithBuildMetadata:
        if configuration.crystal is None:
            raise ValueError("Crystal configuration is required for PointDefectSiteBuilder")

        new_material = MaterialWithBuildMetadata.create(
            {
                "name": configuration.crystal.name,
                "lattice": configuration.crystal.lattice.to_dict(),
                "basis": configuration.crystal.basis.to_dict(),
            }
        )
        elements = configuration.crystal.basis.elements.values
        new_material.basis.remove_atoms_by_elements(elements)
        new_material.basis.add_atom(
            element=configuration.element.chemical_element,
            coordinate=configuration.coordinate,
        )
        return new_material
