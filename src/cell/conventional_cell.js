import { LATTICE_TYPE } from "../lattice/types";
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
    [LATTICE_TYPE.CUB]: unitMatrix,
    [LATTICE_TYPE.FCC]: [
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ],
    [LATTICE_TYPE.BCC]: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ],
    [LATTICE_TYPE.TET]: unitMatrix,
    [LATTICE_TYPE.BCT]: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ],
    [LATTICE_TYPE.ORC]: unitMatrix,
    [LATTICE_TYPE.ORCF]: [
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ],
    [LATTICE_TYPE.ORCI]: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ],
    [LATTICE_TYPE.ORCC]: [
        [1, -1, 0],
        [1, 1, 0],
        [0, 0, 1],
    ],
    [LATTICE_TYPE.HEX]: unitMatrix,
    [LATTICE_TYPE.RHL]: unitMatrix,
    [LATTICE_TYPE.MCL]: unitMatrix,
    [LATTICE_TYPE.MCLC]: [
        [1, -1, 0],
        [1, 1, 0],
        [0, 0, 1],
    ],
    [LATTICE_TYPE.TRI]: unitMatrix,
    [`${LATTICE_TYPE.TRI}alt`]: unitMatrix,
};

export const PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES = {
    // PRIMITIVE    =>  CONVENTIONAL
    [LATTICE_TYPE.CUB]: LATTICE_TYPE.CUB,
    [LATTICE_TYPE.FCC]: LATTICE_TYPE.CUB,
    [LATTICE_TYPE.BCC]: LATTICE_TYPE.CUB,
    [LATTICE_TYPE.TET]: LATTICE_TYPE.TET,
    [LATTICE_TYPE.BCT]: LATTICE_TYPE.TET,
    [LATTICE_TYPE.ORC]: LATTICE_TYPE.ORC,
    [LATTICE_TYPE.ORCF]: LATTICE_TYPE.ORC,
    [LATTICE_TYPE.ORCI]: LATTICE_TYPE.ORC,
    [LATTICE_TYPE.ORCC]: LATTICE_TYPE.ORC,
    [LATTICE_TYPE.HEX]: LATTICE_TYPE.HEX,
    [LATTICE_TYPE.RHL]: LATTICE_TYPE.RHL,
    [LATTICE_TYPE.MCL]: LATTICE_TYPE.MCL,
    [LATTICE_TYPE.MCLC]: LATTICE_TYPE.MCL,
    [LATTICE_TYPE.TRI]: LATTICE_TYPE.TRI,
    // TODO: Legacy `TRI_alt` type, assert not used and remove
    [`${LATTICE_TYPE.TRI}alt`]: LATTICE_TYPE.TRI,
};

export function isConventionalCellSameAsPrimitiveForLatticeType(latticeType) {
    const multiplier = PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS[latticeType || LATTICE_TYPE.TRI];
    return multiplier === unitMatrix;
}

/**
 * Returns lattice vectors of a conventional cell for a primitive lattice.
 * @param {Lattice} lattice - Lattice instance.
 * @return {Array[]} Cell.vectorsAsArray
 */
export function getConventionalCellFromPrimitiveLattice(lattice) {
    const multiplier = PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS[lattice.type || LATTICE_TYPE.TRI];
    return math.multiply(math.matrix(multiplier), math.matrix(lattice.vectorArrays)).toArray();
}
