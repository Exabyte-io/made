from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.point.atom_at_coordinate.configuration import AtomAtCoordinateConfiguration
from mat3ra.made.tools.build.vacuum.builders import VacuumBuilder
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration


class AtomAtCoordinateBuilder(VacuumBuilder):
    _ConfigurationType = AtomAtCoordinateConfiguration

    def _generate(self, configuration: AtomAtCoordinateConfiguration) -> MaterialWithBuildMetadata:
        if not configuration.crystal:
            raise ValueError("Crystal must be provided for AtomAtCoordinateConfiguration")
        vacuum_configuration = VacuumConfiguration(
            crystal=configuration.crystal,
            size=configuration.crystal.lattice.c,
            direction=AxisEnum.z,
        )
        new_material = super()._generate(vacuum_configuration)
        new_material.basis.add_atom(
            element=configuration.element.chemical_element.value,
            coordinate=configuration.coordinate,
        )
        return new_material
