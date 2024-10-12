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
defect_creation_point_vacancy.py

defect_creation_point_substitution.py

defect_creation_point_interstitial.py

defect_creation_point_defect_pair.py

3.2. Surface Defects
defect_creation_surface_adatom.py

defect_creation_surface_island.py

defect_creation_surface_terrace.py

3.3. Planar Defects
defect_creation_planar_grain_boundary_slab.py

defect_creation_planar_grain_boundary_bulk.py

4. Passivation
4.1. Surface Passivation
slab_add_passivation.py

4.2. Edge Passivation
monolayer_nanoribbon_add_passivation.py

5. Perturbations
slab_add_perturbation_wave_longitudinal.py

slab_add_perturbation_wave_radial.py
slab_add_perturbation_custom.py


## 1. Ontology

One can think about the process of a creation of a new materials as a workflow, where there is input material, a set of operations that are applied to it, and the output material(s) - Result(s). Following this, we identify the main concepts as:

- (Input) Materials: the "holders" of the structural information (e.g. Silicon FCC crystal).
- Configurations: describe the physical properties of the Result (e.g. defect, interface, etc.).
- Builders: the workflows, transformations applied to Configurations to generate Result(s) .
- BuilderParameters: specific parameters of the digital representation of the transformations.
- (Result) Materials: the output of the Builders (e.g. point defect in Silicon, graphene nanoribbon etc.).

[//]: # (a diagram that explains the relationships between entities)
![Made Tools](images/made-tools.png)

<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 300" width="800" height="300">
  <defs>
    <style>
      .background { fill: #1e1e1e; }
      .section { fill: #2d2d2d; stroke: #3a3a3a; stroke-width: 2; }
      .box { stroke: #505050; stroke-width: 2; }
      .material { fill: #264f78; }
      .feature { fill: #1e4d2b; }
      .builder { fill: #2d2d2d; }
      .metadata { fill: #3c3c3c; }
      .text { font-family: Arial, sans-serif; font-size: 14px; fill: #d4d4d4; }
      .arrow { fill: none; stroke: #d4d4d4; stroke-width: 2; }
    </style>
  </defs>

  <!-- Background -->
  <rect width="100%" height="100%" class="background"/>

  <!-- Input Section -->
  <rect x="10" y="10" width="250" height="280" class="section"/>
  <text x="135" y="35" class="text" text-anchor="middle">Input</text>
  <rect x="20" y="50" width="230" height="230" rx="20" ry="20" class="box"/>
  <text x="135" y="75" class="text" text-anchor="middle">Configuration</text>
  
  <rect x="30" y="90" width="210" height="30" class="material"/>
  <text x="135" y="110" class="text" text-anchor="middle">Material</text>
  
  <rect x="30" y="130" width="210" height="30" class="material"/>
  <text x="135" y="150" class="text" text-anchor="middle">Material</text>
  
  <text x="135" y="190" class="text" text-anchor="middle">...</text>
  
  <rect x="30" y="220" width="210" height="50" class="feature"/>
  <text x="135" y="250" class="text" text-anchor="middle">Feature Description</text>

  <!-- Workflow Section -->
  <rect x="270" y="10" width="250" height="280" class="section"/>
  <text x="395" y="35" class="text" text-anchor="middle">Workflow</text>
  
  <rect x="280" y="50" width="230" height="60" rx="20" ry="20" class="box"/>
  <text x="395" y="85" class="text" text-anchor="middle">Builder Parameters</text>
  
  <rect x="280" y="160" width="230" height="100" class="box"/>
  <text x="395" y="215" class="text" text-anchor="middle">Builder</text>

  <!-- Output Section -->
  <rect x="530" y="10" width="250" height="280" class="section"/>
  <text x="655" y="35" class="text" text-anchor="middle">Output</text>
  <rect x="540" y="50" width="230" height="230" rx="20" ry="20" class="box"/>
  <text x="655" y="75" class="text" text-anchor="middle">(Result)</text>
  
  <rect x="550" y="90" width="210" height="30" class="material"/>
  <text x="655" y="110" class="text" text-anchor="middle">Material</text>
  
  <rect x="550" y="120" width="210" height="30" class="metadata"/>
  <text x="655" y="140" class="text" text-anchor="middle">Metadata</text>
  
  <rect x="550" y="160" width="210" height="30" class="material"/>
  <text x="655" y="180" class="text" text-anchor="middle">Material</text>
  
  <rect x="550" y="190" width="210" height="30" class="metadata"/>
  <text x="655" y="210" class="text" text-anchor="middle">Metadata</text>
  
  <text x="655" y="250" class="text" text-anchor="middle">...</text>

  <!-- Arrows -->
  <path d="M250 180 H280" class="arrow" marker-end="url(#arrowhead)"/>
  <path d="M395 110 V160" class="arrow" stroke-dasharray="5,5" marker-end="url(#arrowhead)"/>
  <path d="M510 210 H540" class="arrow" marker-end="url(#arrowhead)"/>

  <!-- Arrowhead marker -->
  <defs>
    <marker id="arrowhead" markerWidth="10" markerHeight="7" refX="0" refY="3.5" orient="auto">
      <polygon points="0 0, 10 3.5, 0 7" fill="#d4d4d4"/>
    </marker>
  </defs>
</svg>

for a specific case of the creation of an interface:

[//]: # (a diagram that explains the relationships for an Interface)
![Made Tools](images/made-tools.png)


## 2. Code Structure.

### 2.1. The Approach.

The approach follows the ontology described above. The main classes/concepts are:

- Builders and BuilderParameters.
- Configurations.
- Factories: methods that manage the creation of different types of builders.


### 2.2. Top-level Entrypoint Functions.

For usability purposes, we provide top-level functions that allow for the creation of materials with a single call. These functions are:

[//]: # (Explain the helper functions and their purpose `create_interface`.)

- `create_interface`
- `create_defect`
- `create_perturbation`



## 3. Usage

### 3.1. Installation

See top-level [README](LINK) for installation instructions.

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
print(perturbed_slab.to_json())
```



## Links

1. [Exabyte Source of Schemas and Examples (ESSE), Github Repository](https://github.com/exabyte-io/exabyte-esse)
2. [ECMAScript 2015 Language Specifications](https://www.ecma-international.org/ecma-262/6.0/)
3. [POSCAR file format, official website](https://cms.mpi.univie.ac.at/vasp/guide/node59.html)
4. [XYZ file format, Wikipedia](https://en.wikipedia.org/wiki/XYZ_file_format)
5. [Quantum ESPRESSO, Official Website](https://www.quantum-espresso.org/)
6. [Mocha, Official Website](https://mochajs.org/)
7. [Git LFS, Official Website](https://git-lfs.github.com/)
