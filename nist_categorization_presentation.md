## 1  The Challenge

Most existing datasets for ML training assume materials are **perfect 3‑D bulk crystals**.
This mismatch creates three pain points:

1. **Limited realism** – Most open datasets and codes still assume **perfect 3‑D bulk crystals**, yet surfaces, interfaces, thin films, and defects that dominate device behavior.
2. **Jargon overload** – "monolayer", "slab", "heterostack", etc., mean different things across sub‑fields, slowing collaboration, 
3. **Re‑implementation fatigue** – every lab writes one‑off scripts to build the same structures, hurting reproducibility and maintenance.

---

## 2  Project Goal

Create an **open, extensible ecosystem** that can:

1. **Categorize** realistic material entities (2‑D, 1‑D, 0‑D, defective, processed).
2. **Generate** these structures programmatically for simulation and ML.
3. **Standardize** data exchange across tools, labs, and federated databases.

---

## 3  Why a FAIR Category System?

Our categorization embeds FAIR principles by design:

* **Findable** – Each material instance carries a persistent *Category ID* + hashable parameters.
* **Accessible** – JSON schemas and example files are openly hosted on PyPI & GitHub.
* **Interoperable** – The same schema drives both a Python runtime and machine‑readable REST outputs.
* **Reusable** – Builders re‑generate structures bit‑for‑bit, documenting provenance and parameters.

This aligns directly with NIST's **FAIR Digital Object (FDO) Framework** initiative.

---

## 4  Two‑Sided Approach to Categorization

We introduce a **two‑sided categorization system** that acts as a Rosetta Stone between domain scientists and software engineers:

* **Scientific side** – names and parameters exactly as used in publications (Miller indices, terminations, twist angles, vacancy types).
* **Software side** – object‑oriented classes that embody those same concepts so any codebase can build, validate, and transform structures reproducibly.


| Perspective    | What it captures                     | Example                                  |
| -------------- | ------------------------------------ | ---------------------------------------- |
| **Scientific** | Terms used in papers & lab notebooks | `slab`, `miller_indices`, `vacancy`      |
| **Software**   | Executable Python objects & methods  | `class Slab(...)`, `VacancyBuilder(...)` |

A shared **JSON‑Schema ID** binds both layers, turning any *methods section* into an *actionable digital object*.


This dual view links a paper's methods section directly to an executable JSON description.

---

## 4  Two‑Sided Approach — Details

| Perspective    | Core Question                              | Representation                                               |
| -------------- | ------------------------------------------ | ------------------------------------------------------------ |
| **Scientific** | "How do I cite this structure in a paper?" | Ontology terms (`monolayer`, `miller_indices`, `vacancy`).   |
| **Software**   | "How does my script build or parse it?"    | Python/TS classes (`Monolayer`, `MillerIndices`, `Vacancy`). |

Both sides reference the **same JSON‑Schema IDs**, so updates propagate automatically.

---

## 5  Digital Stack Components

1. **Data Standards** — JSON Schemas & worked examples
   • Core ontology (`mat3ra‑esse`) and a library of canonical example files.
2. **Initial Data** — curated structures & workflow templates
   • Bulk crystals, slabs, interfaces, plus ready‑to‑run Jupyter workflow notebooks.
3. **Runtime Code** — functions driven by the standards
   • Builders, converters, and validators implemented in `mat3ra‑made`.
4. **Applications** — user‑facing layers
   • Notebook gallery (Jupyter‑Lite), Materials Designer web‑app, and CI pipelines that validate every contributed POSCAR.

---


## 6  Interactive Notebooks & Deployment

* **15 + notebooks** showcase typical tasks: slab builder, interface matching, defect insertion.
* Hosted at **jupyterlite.mat3ra.com** — run in‑browser, no install required.
* Integrated into the **Materials Designer** web‑app for point‑and‑click workflows.

---

## 7  Ontology in Practice

Example JSON configuration for a 2D slab structure:

```jsonc
{
  "material_category": "pristine_structures/two_dimensional/slab",
  "builders": [
    {
      "type": "StructureBuilder.Slab",
      "parameters": {
        "crystal": {
          "name": "Silicon",
          "basis": {
            "elements": [{"id": 0, "value": "Si"}, {"id": 1, "value": "Si"}],
            "coordinates": [{"id": 0, "value": [0.0, 0.0, 0.0]}, {"id": 1, "value": [0.25, 0.25, 0.25]}]
          },
          "lattice": {"a": 5.43, "b": 5.43, "c": 5.43, "alpha": 90, "beta": 90, "gamma": 90}
        },
        "miller_indices": [1, 0, 0],
        "number_of_layers": 4,
        "vacuum": 15.0
      }
    },
    {
      "type": "ConfigurationBuilder.AtomicLayers",
      "parameters": {
        "termination": "Si_P4/mmm_1",
        "repetitions": 4
      }
    },
    {
      "type": "ConfigurationBuilder.Vacuum", 
      "parameters": {
        "size": 15.0,
        "direction": "z"
      }
    }
  ]
}
```

---


## 8  Category Tree (condensed)

### 8.1  Entities

* **3‑D** | `crystal`, `void`, `supercell`, `strained_uniform`
* **2‑D** | `slab_unit_cell`, `monolayer`, `atomic_layers_unique`
* **1‑D** | `ribbon_section`, `nanowire`
* **0‑D** | `nanoparticle`, `vacancy`, `interstitial`

### 8.2  Operations

`strain`, `repeat`, `stack`, `merge`, `perturb`, `fold`

### 8.3  Material Groups

`pristine_structures`, `compound_pristine_structures`, `defective_structures`, `processed_structures`

*The full tree with 240 leaf categories is provided in `categories_tree.yaml` in the repository.*

---

## 9 
| Step | Action                                       | Category Tag                             |
| ---- | -------------------------------------------- |------------------------------------------|
| 1    | `Materials.get_by_name_first_match("Ni")`    | `crystal`                                |
| 2    | `SlabConfiguration.from_parameters(...)`     | `slab_unit_cell`                         |
| 3    | `ZSLInterfaceAnalyzer(...)`                  | `interface_analysis/zsl_matching`        |
| 4    | `analyzer.get_strained_configurations()`     | `strained_configurations`                |
| 5    | `SlabStrainedSupercellWithGapConfiguration`  | `interface_gap_configuration`            |
| 6    | `create_interface(interface_config)`         | `compound_pristine_structures/interface` |

**Key Features:**
- **ZSL Algorithm**: Finds optimal supercell matches with minimal strain
- **Gap Control**: Precise interlayer distance control (3.0 Å)
- **Strain Analysis**: Automatic calculation of lattice mismatch and strain matrices
- **Material Categories**: Each step maps to specific ontology categories

Resulting interface: ~5 atoms, optimized lattice matching with < 2% strain.

---

## 10  Collaboration Opportunity with NIST

**How this accelerates the FAIR Digital Object Framework:**

1. **Rich Test Content** – A comprehensive category library converts directly into candidate FDOs, ready for registry pilots.
2. **Reference Implementation** – An open‑source Python runtime shows how FDO metadata drives real‑world builders and converters.
3. **Cross‑walks** – JSON‑LD contexts map our schema to NIST MGI, Materials Project, and other ontologies.
4. **Validation Service** – Mat3ra's CI can serve as a public endpoint to verify FDO compliance of uploaded POSCAR‑like files.

*Next step*: co‑author a short white‑paper and pilot registry demo with Dr. Zachary Trautt's group.

---

## 11  Conclusions

* A dense, machine‑readable category system unites **academia**, **R&D**, and **automation**.
* Open schemas and builders to cut boilerplate.
* Community input is essential — contribute new edge cases!


## 12  Appendix: Example Interface Configuration

** Example JSON configuration for a Graphene-Nickel interface: **

```json
{
  "type": "InterfaceConfiguration",
  "stack_components": [
    {
      "type": "SlabStrainedSupercellWithGapConfiguration",
      "gap": 3.0,
      "strain_matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
      "xy_supercell_matrix": [[1, 0], [0, 1]],
      "stack_components": [
        {
          "type": "AtomicLayersUniqueRepeatedConfiguration", 
          "crystal": {
            "name": "Ni",
            "basis": {"elements": [{"id": 0, "value": "Ni"}], "...": "..."},
            "lattice": {"a": 3.52, "...": "..."}
          },
          "miller_indices": [1, 1, 1],
          "number_of_repetitions": 3
        },
        {
          "type": "VacuumConfiguration",
          "size": 0.0,
          "direction": "z"
        }
      ]
    },
    {
      "type": "SlabStrainedSupercellWithGapConfiguration", 
      "gap": 3.0,
      "strain_matrix": [[0.964, 0.0, 0.0], [0.0, 0.964, 0.0], [0.0, 0.0, 1.0]],
      "xy_supercell_matrix": [[1, 0], [0, 1]],
      "stack_components": [
        {
          "type": "AtomicLayersUniqueRepeatedConfiguration",
          "crystal": {
            "name": "Graphene", 
            "basis": {"elements": [{"id": 0, "value": "C"}, {"id": 1, "value": "C"}], "...": "..."},
            "lattice": {"a": 2.467, "gamma": 120, "...": "..."},
          },
          "miller_indices": [0, 0, 1],
          "number_of_repetitions": 1
        },
        {
          "type": "VacuumConfiguration",
          "size": 0.0, 
          "direction": "z"
        }
      ]
    },
     {
        "type": "VacuumConfiguration",
        "size": 10.0,
        "direction": "z"
     }
  ],
  "direction": "z",
  "xy_shift": [0.0, 0.0]
}
```
## 13 Appendix: Example Interface Building Process

Usage Example — Graphene ∥ Ni(111) Interface with ZSL

An example of creating a Graphene-Nickel interface using ZSL lattice matching:

```python
from mat3ra.standata.materials import Materials
from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from mat3ra.made.tools.analyze.interface.zsl import ZSLInterfaceAnalyzer
from mat3ra.made.tools.build.interface import InterfaceConfiguration, create_interface
from mat3ra.made.tools.build.slab.configurations import SlabStrainedSupercellWithGapConfiguration

# Step 1: Get materials from Standata
graphene = Material.create(Materials.get_by_name_first_match("Graphene"))
nickel = Material.create(Materials.get_by_name_first_match("Ni"))

# Step 2: Create slab configurations  
substrate_slab_config = SlabConfiguration.from_parameters(
    nickel, 
    miller_indices=(1, 1, 1), 
    number_of_layers=3, 
    vacuum=0.0
)

film_slab_config = SlabConfiguration.from_parameters(
    graphene, 
    miller_indices=(0, 0, 1), 
    number_of_layers=1, 
    vacuum=0.0
)

# Step 3: Analyze interface with ZSL to find optimal matching
analyzer = ZSLInterfaceAnalyzer(
    substrate_slab_configuration=substrate_slab_config,
    film_slab_configuration=film_slab_config,
    max_area=50.0,
)

# Step 4: Get strained configurations from ZSL analysis
strained_configs = analyzer.get_strained_configurations()
selected_match = strained_configs[0]  # Choose best match

# Step 5: Add gap between substrate and film
gap = 3.0  # Angstroms
substrate_config_with_gap = SlabStrainedSupercellWithGapConfiguration(
    **selected_match.substrate_configuration.to_dict(), gap=gap
)
film_config_with_gap = SlabStrainedSupercellWithGapConfiguration(
    **selected_match.film_configuration.to_dict(), gap=gap
)

# Step 6: Create interface configuration and build interface
interface_config = InterfaceConfiguration(
    stack_components=[substrate_config_with_gap, film_config_with_gap]
)

interface = create_interface(interface_config)
```
