[![npm version](https://badge.fury.io/js/%40mat3ra%2Fmade.svg)](https://badge.fury.io/js/%40mat3ra%2Fmade)
[![License: Apache](https://img.shields.io/badge/License-Apache-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

# Made

Made is a library for **MA**terials **DE**sign. It allows for creating and manipulating material structures from atoms up. The library is aimed to be used for the development of web applications, both on the client (web browser) and server (eg. Node.js) side. It has implementations in Python (including Pyodide) and JavaScript/TypeScript.

## 1. Overview

The package provides software concepts for interacting with Materials-related data structures. The concepts follow the ESSE Data Convention [[1]](#links).

## 2. Installation

### 2.1. JavaScript/TypeScript

From NPM for use within a software project:

```bash
npm install @mat3ra/made

```

### 2.2. Python

From PyPI for use within a software project:

```bash
pip install mat3ra-made
```

When willing to use the optional `tools` module, install the package with the following command:

```bash
pip install "mat3ra-made[tools]"
```


## 3. Functionality

As below

- High-level classes for the representation of the [Material](src/material.js) and the corresponding structural information, ie:
    - [Basis](src/basis/basis.js),
    - [Lattice](src/lattice/lattice.js),
    - [ReciprocalLattice](src/lattice/reciprocal/lattice_reciprocal.js),
    - [Cell](src/cell/cell.js),
    - [AtomicConstraints](src/constraints/constraints.js)
    - and others to be added.
- input/output support, including:
    - POSCAR [[3]](#links),
    - XYZ [[4]](#links),
    - Quantum ESPRESSO [[5]](#links),
    - and others to be added.
- structural generation and analysis tools:
    - [supercell](src/tools/supercell.js)
    - [surfaces](src/tools/surface.js)
    - [combinatorial sets](src/parsers/xyz_combinatorial_basis.js)
    - [interpolated sets for chemical reactions](src/tools/basis.js)


## 4. Contribution

This repository is an [open-source](LICENSE.md) work-in-progress and we welcome contributions.

We regularly deploy the latest code containing all accepted contributions online as part of the [Mat3ra.com](https://mat3ra.com) platform, so contributors will see their code in action there.

We suggest forking this repository and introducing the adjustments there to be considered for merging into this repository as explained in more details [here](https://gist.github.com/Chaser324/ce0505fbed06b947d962), for example.

### 4.1. Source code conventions

Object-oriented design patterns encapsulate key concepts following the conventions below.

1. Classes follow the Exabyte Data Convention and data structures defined in ESSE [[1]](#links)

2. Only materials-related code is considered. Properties related to [simulation model](https://docs.exabyte.io/models/overview/) parameters (eg. type of approximation, numerical parameters) shall go elsewhere.

3. `tools` directory contains helper functions that act on one or more classes and include an external parameter. Functions that use class data without any external parameters should be implemented inside the class. For example, `basis.clone()` is implemented in `Basis`, but basis repetition is implemented as a tool in the correspondingly named function ([tools/basis.js#repeat](src/tools/basis.js)) because the repetion requires a parameter external to basis - number of repetitions in 3 spatial dimensions.

4. [Deprecated, use mat3ra-parsers or @mat3ra/parsers] `parsers` directory contains the parsers to- and from- ESSE format mentioned in 1. All functionality related to external data conversion is contained in this directory.


### 4.2. TODO list

[Outdated] Desirable features for implementation:

- identify primitive / conventional structures
- support for molecular geometries
- support for polymer geometries
- radial correlation function calculation
- generation of complex atomic shapes:
    - fullerene
    - nanotube
    - nanowire
    - nano-cluster
    - a combination of the above
    - arbitrary atomic arrangement

## 5. Development

### 5.1. JavaScript/TypeScript

#### 5.1.1. Tests

Made tests are written based on Mocha [6](#links) testing framework and can be executed as follows.

```bash
git pull
git lfs pull
```
to get the latest test fixtures from LFS, and then:

```bash
npm install
npm test
```

#### 5.1.2. Important Notes

1. Keep the tests directory structure similar to the main codebase directory structure. Every JS module in the main codebase should have a corresponding module in tests directory which implements the tests for provided functionality.

2. Add tests fixtures into [fixtures](./tests/fixtures) directory. The fixtures are automatically stored on Git LFS [7](#links).

3. If the fixtures are going to be used inside multiple cases, read and export them inside [enums](./tests/enums.js) to avoid code duplicates.

4. [Tests setup module](./tests/setup.js) can be used to implement the hooks that are used to prepare the tests environment.

#### 5.1.3. Using Linter

Linter setup will prevent committing files that don't adhere to the code standard. It will
attempt to fix what it can automatically prior to the commit in order to reduce diff noise. This can lead to "unexpected" behavior where a
file that is staged for commit is not identical to the file that actually gets committed. This happens
in the `lint-staged` directive of the `package.json` file (by using a `husky` pre-commit hook). For example,
if you add extra whitespace to a file, stage it, and try to commit it, you will see the following:

```bash
➜  made git:(feature/SOF-4398-TB) ✗ git add src/basis/constrained_basis.js
➜  made git:(feature/SOF-4398-TB) ✗ git commit -m "Test commit non-linted code"
✔ Preparing...
✔ Running tasks...
✖ Prevented an empty git commit!
✔ Reverting to original state because of errors...
✔ Cleaning up...

  ⚠ lint-staged prevented an empty git commit.
  Use the --allow-empty option to continue, or check your task configuration

husky - pre-commit hook exited with code 1 (error)
```

The staged change may remain but will not have been committed. Then it will look like you still have a staged
change to commit, but the pre-commit hook will not actually commit it for you, quite frustrating! Styling can
be applied manually and fixed by running:

```bash
npm run lint:fix
```

In which case, you may need to then add the linter edits to your staging, which in the example above, puts the
file back to identical with the base branch, resulting in no staged changes whatsoever.

#### 5.1.4. Configuring WebStorm for use with Linter

In order for the WebStorm IDE to take full advantage of the linting configuration, it can be configured in the project:

- `Preferences -> Languages & Frameworks -> JavaScript -> Code Quality Tools -> ESLint`
- Check `Automatic ESLint configuration` which should infer all the configurations from the project directory

### 5.2. Python

#### 5.2.1. Tests

Python 3.8+ is required to run the tests. We recommend using PyEnv to manage Python versions. Tests are written based on PyTest and can be executed as follows.

```bash
virtualenv .venv
source .venv/bin/activate
pip install ".[tests]"
pytest tests/py
```

#### 5.2.2. Important Notes

Conventions:

- The "tools" module has external dependencies on "pymatgen" and "ase" packages and so is meant as optional. When implementing new functionality, the use of ASE is recommended over pymatgen for compatibility purposes.

#### 5.2.3. Developing locally for pyodide

To build and serve locally, use the following command:

```bash
wheel_server
```
More details can be found in the [script documentation](https://github.com/Exabyte-io/utils/blob/main/README.md).

### 5.3. Known Issues

#### 5.3.1. JavaScript/TypeScript

To be added.

#### 5.3.2. Python

As below:

- Python 3.8 tests are failing on Apple Silicon due to https://github.com/materialsproject/pymatgen/issues/3521

## Links

1. [Exabyte Source of Schemas and Examples (ESSE), Github Repository](https://github.com/exabyte-io/exabyte-esse)
2. [ECMAScript 2015 Language Specifications](https://www.ecma-international.org/ecma-262/6.0/)
3. [POSCAR file format, official website](https://cms.mpi.univie.ac.at/vasp/guide/node59.html)
4. [XYZ file format, Wikipedia](https://en.wikipedia.org/wiki/XYZ_file_format)
5. [Quantum ESPRESSO, Official Website](https://www.quantum-espresso.org/)
6. [Mocha, Official Website](https://mochajs.org/)
7. [Git LFS, Official Website](https://git-lfs.github.com/)
