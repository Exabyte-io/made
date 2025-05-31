from typing import List, Union, Optional, Any

from mat3ra.esse.models.material.primitive.two_dimensional.miller_indices import MillerIndicesSchema
from mat3ra.esse.models.material.reusable.two_dimensional.atomic_layers import AtomicLayersSchema
from mat3ra.esse.models.material.reusable.two_dimensional.atomic_layers_unique import (
    AtomicLayersUniqueSchema,
)
from mat3ra.esse.models.material.reusable.two_dimensional.atomic_layers_unique_repeated import (
    AtomicLayersUniqueRepeatedSchema,
)
from mat3ra.esse.models.material.reusable.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import (
    SlabConfigurationSchema,
    AxisEnum,
    VacuumConfigurationSchema,
)

from mat3ra.made.material import Material
from .termination import Termination
from .. import BaseConfigurationPydantic
from ...analyze.lattice import LatticeMaterialAnalyzer
from ...convert import to_pymatgen, from_pymatgen
from ...modify import wrap_to_unit_cell
from ...third_party import PymatgenSlab
from ...third_party import PymatgenSlabGenerator, label_pymatgen_slab_termination


class CrystalLatticePlanes(CrystalLatticePlanesSchema):

    crystal: Material
    _pymatgen_slabs: List[PymatgenSlab]

    @classmethod
    def from_dict(cls, crystal: dict, **kwargs: Any) -> "CrystalLatticePlanes":
        crystal_material = Material(**crystal)
        return cls(
            crystal=crystal_material,
            **kwargs,
        )

    def __init__(
        self, crystal: Material, miller_indices: MillerIndicesSchema, use_conventional_cell=True, **kwargs: Any
    ):
        if use_conventional_cell:
            lattice_analyzer = LatticeMaterialAnalyzer(crystal, use_cartesian=False)
            crystal = lattice_analyzer.material_with_conventional_lattice
        super().__init__(
            crystal=crystal, miller_indices=miller_indices, use_conventional_cell=use_conventional_cell, **kwargs
        )
        self._pymatgen_slabs = self._generate_pymatgen_slabs()

    def _generate_pymatgen_slabs(
        self,
        min_slab_size=1,
        min_vacuum_size=0,
        in_unit_planes: bool = True,
        make_primitive: bool = False,
        symmetrize: bool = False,
    ) -> List[PymatgenSlab]:  # type: ignore
        generator = PymatgenSlabGenerator(
            initial_structure=to_pymatgen(self.crystal),
            miller_index=self.miller_indices,
            min_slab_size=min_slab_size,
            min_vacuum_size=min_vacuum_size,
            in_unit_planes=in_unit_planes,
            primitive=make_primitive,
        )
        raw_slabs = generator.get_slabs(
            # We need to preserve symmetric slabs for different terminations at the surface
            symmetrize=symmetrize
        )

        return raw_slabs

    def get_terminations(self):
        return [Termination.from_string(label_pymatgen_slab_termination(slab)) for slab in self._pymatgen_slabs]

    def get_slabs(
        self,
        min_slab_size: Union[float, int],
        min_vacuum_size: Union[float, int],
        in_unit_planes: bool,
        make_primitive: bool,
        symmetrize: bool,
        use_orthogonal_c: bool,
    ) -> List[Material]:

        pymatgen_slabs = self._generate_pymatgen_slabs(
            min_slab_size=min_slab_size,
            min_vacuum_size=min_vacuum_size,
            in_unit_planes=in_unit_planes,
            make_primitive=make_primitive,
            symmetrize=symmetrize,
        )
        orthogonal_slabs = (
            [slab.get_orthogonal_c_slab() for slab in pymatgen_slabs] if use_orthogonal_c else pymatgen_slabs
        )
        material_slabs = [Material.create(from_pymatgen(slab)) for slab in orthogonal_slabs]
        wrapped_material_slabs = [wrap_to_unit_cell(material) for material in material_slabs]
        return wrapped_material_slabs


class AtomicLayers(AtomicLayersSchema, CrystalLatticePlanes):
    pass


class AtomicLayersUnique(AtomicLayersUniqueSchema, CrystalLatticePlanes):
    pass


class AtomicLayersUniqueRepeated(AtomicLayersUniqueRepeatedSchema, CrystalLatticePlanes):
    crystal: Material
    termination_top: Termination


class VacuumConfiguration(VacuumConfigurationSchema):
    """
    Configuration for adding vacuum to a material.

    Args:
        direction (AxisEnum): The direction to add vacuum along (x, y, or z).
        size (float): Size of the vacuum gap in angstroms. Defaults to 10.0.
    """

    type: str = "VacuumConfiguration"
    direction: AxisEnum = VacuumConfigurationSchema.model_fields["direction"].default
    size: float = VacuumConfigurationSchema.model_fields["size"].default


class SlabConfiguration(SlabConfigurationSchema, BaseConfigurationPydantic):
    type: str = "SlabConfiguration"
    stack_components: List[Union[AtomicLayersUniqueRepeated, VacuumConfiguration]]
    xy_supercell_matrix: Optional[List[List[int]]] = SlabConfigurationSchema.model_fields["xy_supercell_matrix"].default
    direction: AxisEnum = AxisEnum.z

    def __init__(self, stack_components: List[Any], **kwargs: Any):
        stack_components_as_class_instances = []
        for component in stack_components:
            if isinstance(component, dict):
                if component.get("number_of_repetitions", None) is not None:
                    stack_components_as_class_instances.append(AtomicLayersUniqueRepeated.from_dict(**component))
                elif component.get("lattice", None) is not None:
                    stack_components_as_class_instances.append(Material(**component))
            else:
                stack_components_as_class_instances.append(component)

        super().__init__(stack_components=stack_components_as_class_instances, **kwargs)

    @property
    def atomic_layers(self) -> Union[AtomicLayersUniqueRepeated, AtomicLayersUnique]:
        return self.stack_components[0]

    @property
    def vacuum(self) -> float:
        if len(self.stack_components) < 2:
            return 0
        return self.stack_components[1].size

    @property
    def bulk(self) -> Material:
        atomic_layers = self.atomic_layers
        return atomic_layers.crystal

    @property
    def miller_indices(self) -> MillerIndicesSchema:
        atomic_layers = self.atomic_layers
        return atomic_layers.miller_indices

    @property
    def number_of_layers(self) -> int:
        atomic_layers = self.atomic_layers
        return atomic_layers.number_of_repetitions

    def get_terminations(self) -> List[Termination]:
        return self.atomic_layers.get_terminations()

    def get_termination_by_index(self, index=0) -> Termination:
        all_terminations = self.atomic_layers.get_terminations()
        if index < 0 or index >= len(all_terminations):
            raise IndexError(f"Termination index {index} is out of range.")
        return all_terminations[index]

    @property
    def termination(self) -> Termination:
        return self.atomic_layers.termination_top

    @property
    def termination_index(self) -> int:
        all_terminations = self.atomic_layers.get_terminations()
        return all_terminations.index(self.termination)

    @property
    def parameters(self):
        return {
            "bulk": self.bulk,
            "miller_indices": self.miller_indices,
            "number_of_layers": self.number_of_layers,
            "vacuum": self.vacuum,
            "xy_supercell_matrix": self.xy_supercell_matrix,
            "use_conventional_cell": self.atomic_layers.use_conventional_cell,
            "termination_index": self.termination_index,
        }

    @classmethod
    def from_parameters(
        cls,
        bulk: Material,
        miller_indices: MillerIndicesSchema = (0, 0, 1),
        number_of_layers: int = 1,
        vacuum: float = 10.0,
        xy_supercell_matrix: List[List[int]] = None,
        use_conventional_cell: bool = True,
        termination_index: int = 0,
        **kwargs,
    ) -> "SlabConfiguration":
        if xy_supercell_matrix is None:
            xy_supercell_matrix = [[1, 0], [0, 1]]

        temp_planes = CrystalLatticePlanes(
            crystal=bulk,
            miller_indices=miller_indices,
            use_conventional_cell=use_conventional_cell,
        )

        available_terminations = temp_planes.get_terminations()
        termination = available_terminations[termination_index]

        atomic_layers = AtomicLayersUniqueRepeated(
            crystal=bulk,
            miller_indices=miller_indices,
            number_of_repetitions=number_of_layers,
            use_conventional_cell=use_conventional_cell,
            termination_top=termination,
        )

        vacuum_config = VacuumConfiguration(
            direction=AxisEnum.z,
            size=vacuum,
        )

        direction = kwargs.pop("direction", AxisEnum.z)

        return cls(
            stack_components=[atomic_layers, vacuum_config],
            xy_supercell_matrix=xy_supercell_matrix,
            direction=direction,
            **kwargs,
        )
