import { LATTICE_TYPE, LatticeType } from "../lattice/types";
import math from "../math";

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
    [LatticeType.CUB]: unitMatrix,
    [LatticeType.FCC]: [
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ],
    [LatticeType.BCC]: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ],
    [LatticeType.TET]: unitMatrix,
    [LatticeType.BCT]: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ],
    [LatticeType.ORC]: unitMatrix,
    [LatticeType.ORCF]: [
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ],
    [LatticeType.ORCI]: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ],
    [LatticeType.ORCC]: [
        [1, -1, 0],
        [1, 1, 0],
        [0, 0, 1],
    ],
    [LatticeType.HEX]: unitMatrix,
    [LatticeType.RHL]: unitMatrix,
    [LatticeType.MCL]: unitMatrix,
    [LatticeType.MCLC]: [
        [1, -1, 0],
        [1, 1, 0],
        [0, 0, 1],
    ],
    [LatticeType.TRI]: unitMatrix,
    [`${LatticeType.TRI}alt`]: unitMatrix,
};

export const PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES = {
    // PRIMITIVE    =>  CONVENTIONAL
    [LatticeType.CUB]: LatticeType.CUB,
    [LatticeType.FCC]: LatticeType.CUB,
    [LatticeType.BCC]: LatticeType.CUB,
    [LatticeType.TET]: LatticeType.TET,
    [LatticeType.BCT]: LatticeType.TET,
    [LatticeType.ORC]: LatticeType.ORC,
    [LatticeType.ORCF]: LatticeType.ORC,
    [LatticeType.ORCI]: LatticeType.ORC,
    [LatticeType.ORCC]: LatticeType.ORC,
    [LatticeType.HEX]: LatticeType.HEX,
    [LatticeType.RHL]: LatticeType.RHL,
    [LatticeType.MCL]: LatticeType.MCL,
    [LatticeType.MCLC]: LatticeType.MCL,
    [LatticeType.TRI]: LatticeType.TRI,
    // TODO: Legacy `TRI_alt` type, assert not used and remove
    [`${LatticeType.TRI}alt`]: LatticeType.TRI,
};

export function isConventionalCellSameAsPrimitiveForLatticeType(latticeType: LatticeType): boolean {
    const multiplier = PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS[latticeType || LatticeType.TRI];
    return multiplier === unitMatrix;
}

