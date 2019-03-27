import math from "../math";
import {LATTICE_TYPE} from "../lattice/types";

/**
 * Routines for calculating conventional cell vectors from primitive cell Bravais parameters.
 * Following Setyawan, W., & Curtarolo, S. (2010). doi:10.1016/j.commatsci.2010.05.010
 */

const unitMatrix = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
];

const CONVENTIONAL_TO_PRIMITIVE_CELL_MULTIPLIERS = {

    [LATTICE_TYPE.CUB]: unitMatrix,

    [LATTICE_TYPE.FCC]: [
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1]
    ],

    [LATTICE_TYPE.BCC]: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0]
    ],

    [LATTICE_TYPE.TET]: unitMatrix,

    [LATTICE_TYPE.BCT]: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0]
    ],

    [LATTICE_TYPE.ORC]: unitMatrix,

    [LATTICE_TYPE.ORCF]: [
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1]
    ],

    [LATTICE_TYPE.ORCI]: [
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0]
    ],

    [LATTICE_TYPE.ORCC]: [
        [1, -1, 0],
        [1, 1, 0],
        [0, 0, 1]
    ],

    [LATTICE_TYPE.HEX]: unitMatrix,

    [LATTICE_TYPE.RHL]: unitMatrix,

    [LATTICE_TYPE.MCL]: unitMatrix,

    [LATTICE_TYPE.MCLC]: [
        [1, -1, 0],
        [1, 1, 0],
        [0, 0, 1]
    ],

    [LATTICE_TYPE.TRI]: unitMatrix,

    [LATTICE_TYPE.TRI + 'alt']: unitMatrix,

};

/**
 * Returns lattice vectors for a conventional cell for a lattice.
 * @param {Lattice} lattice - Lattice instance.
 * @return {Array[]} Cell.vectorsAsArray
 */
export function conventionalCell(lattice) {
    const multiplier = CONVENTIONAL_TO_PRIMITIVE_CELL_MULTIPLIERS[lattice.type || LATTICE_TYPE.TRI];
    return math.multiply(math.matrix(multiplier), math.matrix(lattice.vectorArrays)).toArray();
}
