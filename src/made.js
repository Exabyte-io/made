import {tolerance, coefficients, units, ATOMIC_COORD_UNITS} from "./constants";
import MadeMath from "./math";
import {Material, defaultMaterialConfig} from "./material";
import {Lattice} from "./lattice/lattice";
import {Basis} from "./basis/basis";
import {ReciprocalLattice} from "./lattice/reciprocal/lattice_reciprocal";
import parsers from "./parsers/parsers";
import tools from "./tools/index";
import {DEFAULT_LATTICE_UNITS, LATTICE_TYPE_CONFIGS} from "./lattice/types";
import {ELEMENTS, ELEMENTS_BY_SYMBOL} from "./periodic_table/elements";
import {ArrayWithIds} from "./abstract/array_with_ids";
import {AtomicConstraints} from "./constraints/constraints";

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
    ELEMENTS,
    ELEMENTS_BY_SYMBOL,
    LATTICE_TYPE_CONFIGS,
    DEFAULT_LATTICE_UNITS,
    primitive: {
        ArrayWithIds,
    },
};
