import {tolerance, coefficients, units, atomicCoordUnits} from "./constants";
import MadeMath from "./math";
import {Material, MaterialMixin, defaultMaterialConfig} from "./material";
import {Lattice} from "./lattice/lattice";
import {Basis} from "./basis/basis";
import {ReciprocalLattice} from "./lattice/reciprocal/lattice_reciprocal";
import parsers from "./parsers/parsers";
import tools from "./tools/index";
import {DEFAULT_LATTICE_UNITS, LATTICE_TYPE_CONFIGS} from "./lattice/types";
import {ELEMENTS, ELEMENTS_BY_SYMBOL} from "./periodic_table/elements";
import {ArrayWithIds} from "./primitive";
import {AtomicConstraints} from "./other/constraints";

export const Made = {
    coefficients: coefficients,
    tolerance: tolerance,
    units: units,
    atomicCoordUnits: atomicCoordUnits,
    math: MadeMath,
    Material,
    MaterialMixin,
    defaultMaterial: defaultMaterialConfig,
    Lattice,
    ReciprocalLattice,
    Basis,
    parsers,
    tools,
    ELEMENTS,
    ELEMENTS_BY_SYMBOL,
    LATTICE_TYPE_CONFIGS,
    DEFAULT_LATTICE_UNITS,
    primitive: {
        ArrayWithIds,
    },
    other: {
        AtomicConstraints,
    }
};
