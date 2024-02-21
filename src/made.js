import { ArrayWithIds } from "./abstract/array_with_ids";
import { Basis } from "./basis/basis";
import { ATOMIC_COORD_UNITS, coefficients, tolerance, units } from "./constants";
import { AtomicConstraints } from "./constraints/constraints";
import { Lattice, nonPeriodicLatticeScalingFactor } from "./lattice/lattice";
import { ReciprocalLattice } from "./lattice/reciprocal/lattice_reciprocal";
import { DEFAULT_LATTICE_UNITS, LATTICE_TYPE_CONFIGS } from "./lattice/types";
import { defaultMaterialConfig, Material, MaterialMixin } from "./material";
import MadeMath from "./math";
import parsers from "./parsers/parsers";
import tools from "./tools/index";

export const Made = {
    coefficients,
    tolerance,
    units,
    ATOMIC_COORD_UNITS,
    math: MadeMath,

    Material,
    MaterialMixin,
    defaultMaterialConfig,
    Lattice,
    nonPeriodicLatticeScalingFactor,
    ReciprocalLattice,
    Basis,
    AtomicConstraints,

    parsers,
    tools,
    LATTICE_TYPE_CONFIGS,
    DEFAULT_LATTICE_UNITS,
    primitive: {
        ArrayWithIds,
    },
};
