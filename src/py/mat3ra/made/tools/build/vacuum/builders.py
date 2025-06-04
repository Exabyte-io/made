from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder
from .configuration import VacuumConfiguration


def create_vacuum_material(reference: Material, vacuum: "VacuumConfiguration", direction: str) -> Material:
    # TODO: update to handle other directions
    a_vector, b_vector = reference.lattice.vector_arrays[:2]
    vacuum_lattice = reference.lattice.from_vectors_array(
        [a_vector, b_vector, [0, 0, vacuum.size]], reference.lattice.units, reference.lattice.type
    )
    return Material.create(
        {
            "name": "Vacuum",
            "lattice": vacuum_lattice.to_dict(),
            "basis": {"elements": [], "coordinates": []},
        }
    )


class VacuumBuilder(BaseBuilder):
    _ConfigurationType = VacuumConfiguration

    def get_material(self, configuration: VacuumConfiguration) -> Material:
        reference = configuration.crystal
        size = configuration.size
        direction = configuration.direction
        return create_vacuum_material(reference, self.crystal.size)
