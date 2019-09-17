import {tolerance, coefficients, units, ATOMIC_COORD_UNITS} from "./constants";
import MadeMath from "./math";
import {Material, defaultMaterialConfig} from "./material";
import {Lattice} from "./lattice/lattice";
import {Basis} from "./basis/basis";
import {ReciprocalLattice} from "./lattice/reciprocal/lattice_reciprocal";
import parsers from "./parsers/parsers";
import tools from "./tools/index";
import {DEFAULT_LATTICE_UNITS, LATTICE_TYPE_CONFIGS, LATTICE_TYPE} from "./lattice/types";
import {ArrayWithIds} from "./abstract/array_with_ids";
import {AtomicConstraints} from "./constraints/constraints";

import {SlabModelBuilder} from "./slab/SlabModelBuilder";

export const Made = {
    coefficients,
    tolerance,
    units,
    ATOMIC_COORD_UNITS,
    math: MadeMath,

    Material,
    defaultMaterialConfig,
    Lattice,
    ReciprocalLattice,
    Basis,
    AtomicConstraints,

    parsers,
    tools,
    LATTICE_TYPE,
    LATTICE_TYPE_CONFIGS,
    DEFAULT_LATTICE_UNITS,
    primitive: {
        ArrayWithIds,
    },
    SlabModelBuilder,
};
