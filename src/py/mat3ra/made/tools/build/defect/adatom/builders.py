# TODO: only exact cooridante and equidistant use this, for Crystal site we need SlabStackBuilder
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.other import get_closest_site_id_from_coordinate
from mat3ra.made.tools.build.defect.adatom.configuration import CrystalSiteAdatomConfiguration
from mat3ra.made.tools.build.defect.adatom.configuration import AdatomDefectConfiguration

from mat3ra.made.tools.build.defect.point.builders import PointDefectBuilder
from mat3ra.made.tools.build.defect.slab.builders import SlabStackBuilder
from mat3ra.made.tools.build.defect.slab.helpers import recreate_slab_with_fractional_layers
from mat3ra.made.tools.build.metadata import get_slab_build_configuration
from mat3ra.made.tools.build.slab.helpers import create_slab


class AdatomDefectBuilder(PointDefectBuilder):
    _ConfigurationType = AdatomDefectConfiguration


class CrystalSiteAdatomSlabDefectBuilder(SlabStackBuilder):
    """
    Builder for creating adatom defects using crystal site placement with SlabStackBuilder.

    This builder implements the new approach where:
    1. Original slab (without vacuum)
    2. Added component (single atom in lattice from recreate_slab_with_fractional_layers)
    3. Vacuum layer
    """

    _ConfigurationType = CrystalSiteAdatomConfiguration

    def _generate(self, configuration: CrystalSiteAdatomConfiguration) -> Material:
        # Get the original slab material from the first stack component
        original_slab_material = self._configuration_to_material(configuration.stack_components[0])

        # For now, we'll use a simpler approach: find the closest site to the given coordinate
        # in a slab with additional layers and use that as the resolved coordinate

        # Create a slab with additional layers to provide new crystal sites

        # Get the slab configuration to create additional layers
        slab_configuration = get_slab_build_configuration(original_slab_material.metadata)
        slab_with_additional_layer = create_slab(
            crystal=original_slab_material,
            miller_indices=slab_configuration.atomic_layers.miller_indices,
            termination=slab_configuration.atomic_layers.termination_top,
            number_of_layers=slab_configuration.number_of_layers + 1,
            vacuum=slab_configuration.vacuum,
            xy_supercell_matrix=[[1, 0], [0, 1]],  # Use default supercell matrix
        )

        # Find the closest site in the new slab to the target coordinate
        closest_site_id = get_closest_site_id_from_coordinate(slab_with_additional_layer, configuration.coordinate)
        resolved_coordinate = slab_with_additional_layer.coordinates_array[closest_site_id]

        # Create the added component (single atom layer) using recreate_slab_with_fractional_layers
        added_layer = recreate_slab_with_fractional_layers(original_slab_material, number_of_layers=1)

        # Create a single atom material at the resolved coordinate with the specified element
        single_atom_material = added_layer.clone()
        # Clear all atoms and add just the single atom at the resolved coordinate
        elements = added_layer.basis.elements.values
        single_atom_material.basis.remove_atoms_by_elements(elements)
        single_atom_material.basis.add_atom(element=configuration.element, coordinate=resolved_coordinate)

        # Update the configuration with the new added component
        updated_configuration = configuration.model_copy()
        updated_configuration.stack_components[1] = single_atom_material

        # Use parent class to generate the final stacked material
        return super()._generate(updated_configuration)
