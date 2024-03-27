import { LatticeTypeSchema } from "@mat3ra/esse/dist/js/types";

/**
 * Routines for calculating conventional cell vectors from primitive cell Bravais parameters.
 * Following Setyawan, W., & Curtarolo, S. (2010). doi:10.1016/j.commatsci.2010.05.010
 */

const unitMatrix = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
];

// (Conventional cellVectors) = (Primitive cellVectors) * (PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIER matrix)
export const PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS = {
    // PRIMITIVE    =>  CONVENTIONAL
    CUB: unitMatrix,
    FCC: [
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ],
    BCC: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ],
    TET: unitMatrix,
    BCT: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ],
    ORC: unitMatrix,
    ORCF: [
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ],
    ORCI: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ],
    ORCC: [
        [1, -1, 0],
        [1, 1, 0],
        [0, 0, 1],
    ],
    HEX: unitMatrix,
    RHL: unitMatrix,
    MCL: unitMatrix,
    MCLC: [
        [1, -1, 0],
        [1, 1, 0],
        [0, 0, 1],
    ],
    TRI: unitMatrix,
    TRIalt: unitMatrix,
};

export const PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES: {
    [key in LatticeTypeSchema | "TRIalt"]: LatticeTypeSchema;
} = {
    // PRIMITIVE    =>  CONVENTIONAL
    CUB: "CUB",
    FCC: "CUB",
    BCC: "CUB",
    TET: "TET",
    BCT: "TET",
    ORC: "ORC",
    ORCF: "ORC",
    ORCI: "ORC",
    ORCC: "ORC",
    HEX: "HEX",
    RHL: "RHL",
    MCL: "MCL",
    MCLC: "MCL",
    TRI: "TRI",
    // TODO: Legacy `TRI_alt` type, assert not used and remove
    TRIalt: "TRI",
};

export function isConventionalCellSameAsPrimitiveForLatticeType(
    latticeType: LatticeTypeSchema,
): boolean {
    const multiplier = PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS[latticeType || "TRI"];
    return multiplier === unitMatrix;
}
