[![npm version](https://badge.fury.io/js/%40mat3ra%2Fmade.svg)](https://badge.fury.io/js/%40mat3ra%2Fmade)
[![License: Apache](https://img.shields.io/badge/License-Apache-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

# Tools

Tools for creating and manipulating materials structures.

## 0. Overview

Module provides classes and functions to build and manipulate materials.

The following material types can be generated.

### 0.1. Single-Material Structures

- 3D Structures
  - Slabs
- 2D Structures
  - Monolayers (equivalent to Slabs)
  - Nanoribbons
- (Quasi) Non-periodic Structures
  - Nanostructures (spherical, cylindrical, etc.)

### 0.2. Multi-Material Structures

- Interfaces
  - with no strain matching
  - with strain matching (ZSL)
  - twisted bilayers, with commensurate lattice search
- Heterostructures
  - with no strain matching
  - with strain matching (ZSL)

- Stacked Nanoribbons

### 0.3. Defects

- Point Defects
  - vacancy
  - substitution
  - interstitial
  - defect pair
- Surface Defects
  - adatom
  - island
  - terrace
- Planar Defects
  - grain boundary in a slab
  - grain boundary in a bulk

### 0.4. Passivation

- Surface Passivation
  - Add passivation to slab
- Edge Passivation
  - Add passivation to monolayer nanoribbon

### 0.5. Perturbations

- Add longitudinal wave perturbation to slab
- Add radial wave perturbation to slab
- Add custom perturbation to slab

## 1. Ontology

One can think about the process of a creation of a new materials as a workflow, where there is input material, a set of operations that are applied to it, and the output material(s) - Result(s). Following this, we identify the main concepts as:

- (Input) Materials: the "holders" of the structural information (e.g. Silicon FCC crystal).
- Configurations: describe the physical properties of the Result (e.g. defect, interface, etc.).
- Builders: the workflows, transformations applied to Configurations to generate Result(s) .
- BuilderParameters: specific parameters of the digital representation of the transformations.
- (Result) Materials: the output of the Builders (e.g. point defect in Silicon, graphene nanoribbon etc.).

![Diagram](../../../../../images/made/tools/tools_diagram_general.png)

for a specific case of the creation of an interface:

[//]: # (a diagram that explains the relationships for an Interface)
![Made Tools](../../../../../images/made/tools/tools_diagram_interface.png)


## 2. Code Structure.

### 2.1. The Approach.

The approach follows the ontology described above. The main classes/concepts are:

- Builders and BuilderParameters.
- Configurations.
- Factories: methods that manage the creation of different types of builders.


### 2.2. Top-level Entrypoint Functions.

For usability purposes, we provide top-level functions that allow for the creation of materials with a single call. These functions are:

[//]: # (Explain the helper functions and their purpose `create_interface`.)

The function with plural name takes in Configuration and Builder as arguments and returns the Result Materials list when multiple results are possible.

- `create_interfaces`
- `create_defects`


The function with singular name takes in Configuration and Builder as arguments and returns single Result Material object when only one result is possible, or the first, most optimal one, from the list of results ordered by some criteria.

- `create_interface`
- `create_defect`
- `create_slab_defect`
- `create_slab`
- `create_perturbation`
- `create_passivation`
- `create_nanoribbon`
- `create_grain_boundary`

Some of the top-level functions are helpers that allow to get additional properties from configuration to refine the build process, such as:
- `get_terminations` for slab
- `get_unique_coordination_numbers` for passivation
- `get_optimal_film_displacement` for interface


## 3. Usage

### 3.1. Installation

See top-level [README](../../../../../README.md) for installation instructions.

### 3.2. Usage in Python

Here is an example workflow using made.tools.build:

Define Configuration:

```python
from mat3ra.made.tools.build.perturbation.configuration import PerturbationConfiguration
from mat3ra.made.tools.build.perturbation.builders import SlabPerturbationBuilder

perturbation_config = PerturbationConfiguration(
    material=initial_material,
    perturbation_function_holder=perturbation_function_holder,
    use_cartesian_coordinates=False
)
```

then, create Builder and generate Material:

```python
builder = SlabPerturbationBuilder()
perturbed_slab = builder.get_material(perturbation_config)
```

or, alternatively, use the top-level function:

```python
from mat3ra.made.tools import create_perturbation

perturbed_slab = create_perturbation(perturbation_config,
                                     builder=SlabPerturbationBuilder)
```

Finally, one can review the content of the material:

```python
print(perturbed_slab.to_dict())
```



## Links

1. [Exabyte Source of Schemas and Examples (ESSE), Github Repository](https://arxiv.org/abs/1902.10838, "Exabyte Source of Schemas and Examples (ESSE), Github Repository")
2. [ASE](https://wiki.fysik.dtu.dk/ase/, "ASE")
3. [Pymatgen](https://pymatgen.org/, "Pymatgen")
