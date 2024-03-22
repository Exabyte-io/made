[![npm version](https://badge.fury.io/js/%40mat3ra%2Fmade.svg)](https://badge.fury.io/js/%40mat3ra%2Fmade)
[![License: Apache](https://img.shields.io/badge/License-Apache-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

# Made

Made is a library for **MA**terials **DE**sign in JavaScript/TypeScript and Python. It allows for the creation and manipulation of material structures from atoms up on the web. The library is aimed to be used for the development of web applications, both on the client (web browser) and server (eg. Node.js) side.

The library was originally designed as part of and presently powers materials design capabilities of the [Exabyte.io](https://exabyte.io) platform. For example, [this page](https://platform.exabyte.io/demo/materials/n3HSzCmyoctgJFGGE) representing a crystal of Silicon online uses Made.js.

Exabyte.io believe in a collaborative future of materials design on the web.

## Functionality

As below:

- the package provides a software environment for interacting with Materials-related data structures from ESSE Data Convention [[1]](#links) and is written in ECMAScript 2015 for use on the web
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

The package is written in a modular way easy to extend. Contributions can be in the form of additional tools or modules you develop, or feature requests and [bug/issue reports](https://help.github.com/articles/creating-an-issue/).

## Installation

From NPM for use within a software project:

```bash
npm install @mat3ra/made

```

From source to contribute to development:

```bash
git clone git@github.com:Exabyte-io/made
```

## Contribution

This repository is an [open-source](LICENSE.md) work-in-progress and we welcome contributions.

### Why contribute?

We regularly deploy the latest code containing all accepted contributions online as part of the [Mat3ra.com](https://mat3ra.com) platform, so contributors will see their code in action there.

### Adding new functionality

We suggest forking this repository and introducing the adjustments there to be considered for merging into this repository as explained in more details [here](https://gist.github.com/Chaser324/ce0505fbed06b947d962), for example.

### Source code conventions

Made.js is written in EcmaScript 6th edition [[2]](#links) with the application of object-oriented design patterns encapsulating key concepts following the conventions below.

1. Classes follow the Exabyte Data Convention and data structures defined in ESSE [[1]](#links)

2. Only materials-related code is considered. Properties related to [simulation model](https://docs.exabyte.io/models/overview/) parameters (eg. type of approximation, numerical parameters) shall go elsewhere.

3. `tools` directory contains helper functions that act on one or more classes and include an external parameter. Functions that use class data without any external parameters should be implemented inside the class. For example, `basis.clone()` is implemented in `Basis`, but basis repetition is implemented as a tool in the correspondingly named function ([tools/basis.js#repeat](src/tools/basis.js)) because the repetion requires a parameter external to basis - number of repetitions in 3 spatial dimensions.

4. `parsers` directory contains the parsers to- and from- ESSE format mentioned in 1. All functionality related to external data conversion is contained in this directory. 


### TODO list

Desirable features for implementation:

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

## Tests

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

### Tests Important Notes

1. Keep the tests directory structure similar to the main codebase directory structure. Every JS module in the main codebase should have a corresponding module in tests directory which implements the tests for provided functionality.

2. Add tests fixtures into [fixtures](./tests/fixtures) directory. The fixtures are automatically stored on Git LFS [7](#links).

3. If the fixtures are going to be used inside multiple cases, read and export them inside [enums](./tests/enums.js) to avoid code duplicates.

4. [Tests setup module](./tests/setup.js) can be used to implement the hooks that are used to prepare the tests environment.

## Using Linter

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

## Configuring WebStorm for use with Linter

In order for the WebStorm IDE to take full advantage of the linting configuration, it can be configured in the project:

- `Preferences -> Languages & Frameworks -> JavaScript -> Code Quality Tools -> ESLint`
- Check `Automatic ESLint configuration` which should infer all the configurations from the project directory

## Links

1. [Exabyte Source of Schemas and Examples (ESSE), Github Repository](https://github.com/exabyte-io/exabyte-esse)
2. [ECMAScript 2015 Language Specifications](https://www.ecma-international.org/ecma-262/6.0/)
3. [POSCAR file format, official website](https://cms.mpi.univie.ac.at/vasp/guide/node59.html)
4. [XYZ file format, Wikipedia](https://en.wikipedia.org/wiki/XYZ_file_format)
5. [Quantum ESPRESSO, Official Website](https://www.quantum-espresso.org/)
6. [Mocha, Official Website](https://mochajs.org/)
7. [Git LFS, Official Website](https://git-lfs.github.com/)
