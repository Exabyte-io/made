from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from .configuration import VacuumConfiguration
from ..utils import AXIS_TO_INDEX_MAP


class VacuumBuilder(BaseBuilder):
    _ConfigurationType = VacuumConfiguration

    def get_material(self, configuration: VacuumConfiguration) -> Material:
        reference = configuration.crystal
        size = configuration.size
        direction: AxisEnum = configuration.direction

        lattice_vectors = reference.lattice.vector_arrays.copy()

        # Replace the target direction with the vacuum vector
        axis = AXIS_TO_INDEX_MAP[direction.value]
        vacuum_vector = [0.0, 0.0, 0.0]
        vacuum_vector[axis] = size
        lattice_vectors[axis] = vacuum_vector

        vacuum_lattice = reference.lattice.from_vectors_array(
            lattice_vectors, reference.lattice.units, reference.lattice.type
        )

        return Material.create(
            {
                "name": "Vacuum",
                "lattice": vacuum_lattice.to_dict(),
                "basis": {"elements": [], "coordinates": []},
            }
        )
