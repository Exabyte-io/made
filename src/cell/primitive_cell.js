/* eslint no-unused-vars: 0 */
import { LATTICE_TYPE } from "../lattice/types";
import math from "../math";

/**
 * Routines for calculating primitive cell vectors from conventional cell Bravais parameters.
 * Following Setyawan, W., & Curtarolo, S. (2010). doi:10.1016/j.commatsci.2010.05.010
 */
const PRIMITIVE_CELLS = {
    [LATTICE_TYPE.CUB]: ({ a, b, c, alpha, beta, gamma }) => {
        return [
            [a, 0, 0],
            [0, a, 0],
            [0, 0, a],
        ];
    },

    [LATTICE_TYPE.FCC]: ({ a, b, c, alpha, beta, gamma }) => {
        return [
            [0.0, a / 2, a / 2],
            [a / 2, 0.0, a / 2],
            [a / 2, a / 2, 0.0],
        ];
    },

    [LATTICE_TYPE.BCC]: ({ a, b, c, alpha, beta, gamma }) => {
        return [
            [-a / 2, a / 2, a / 2],
            [a / 2, -a / 2, a / 2],
            [a / 2, a / 2, -a / 2],
        ];
    },

    [LATTICE_TYPE.TET]: ({ a, b, c, alpha, beta, gamma }) => {
        return [
            [a, 0, 0],
            [0, a, 0],
            [0, 0, c],
        ];
    },

    [LATTICE_TYPE.BCT]: ({ a, b, c, alpha, beta, gamma }) => {
        return [
            [-a / 2, a / 2, c / 2],
            [a / 2, -a / 2, c / 2],
            [a / 2, a / 2, -c / 2],
        ];
    },

    [LATTICE_TYPE.ORC]: ({ a, b, c, alpha, beta, gamma }) => {
        return [
            [a, 0, 0],
            [0, b, 0],
            [0, 0, c],
        ];
    },

    [LATTICE_TYPE.ORCF]: ({ a, b, c, alpha, beta, gamma }) => {
        return [
            [0, b / 2, c / 2],
            [a / 2, 0, c / 2],
            [a / 2, b / 2, 0],
        ];
    },

    [LATTICE_TYPE.ORCI]: ({ a, b, c, alpha, beta, gamma }) => {
        return [
            [-a / 2, b / 2, c / 2],
            [a / 2, -b / 2, c / 2],
            [a / 2, b / 2, -c / 2],
        ];
    },

    [LATTICE_TYPE.ORCC]: ({ a, b, c, alpha, beta, gamma }) => {
        return [
            [a / 2, b / 2, 0],
            [-a / 2, b / 2, 0],
            [0, 0, c],
        ];
    },

    [LATTICE_TYPE.HEX]: ({ a, b, c, alpha, beta, gamma }) => {
        return [
            [a / 2, (-a * math.sqrt(3)) / 2, 0],
            [a / 2, (a * math.sqrt(3)) / 2, 0],
            [0, 0, c],
        ];
    },

    [LATTICE_TYPE.RHL]: ({ a, b, c, alpha, beta, gamma }) => {
        const cosAlpha = math.cos((alpha / 180) * math.PI);
        const cosHalfAlpha = math.sqrt((1 / 2) * (1 + cosAlpha));
        const sinHalfAlpha = math.sqrt((1 / 2) * (1 - cosAlpha));
        return [
            [a * cosHalfAlpha, -a * sinHalfAlpha, 0.0],
            [a * cosHalfAlpha, a * sinHalfAlpha, 0.0],
            [
                (a * cosAlpha) / cosHalfAlpha,
                0.0,
                a * math.sqrt(1 - (cosAlpha * cosAlpha) / (cosHalfAlpha * cosHalfAlpha)),
            ],
        ];
    },

    [LATTICE_TYPE.MCL]: ({ a, b, c, alpha, beta, gamma }) => {
        const cosAlpha = math.cos((alpha / 180) * math.PI);
        return [
            [a, 0, 0],
            [0, b, 0],
            [0, c * cosAlpha, c * math.sqrt(1 - cosAlpha * cosAlpha)],
        ];
    },

    [LATTICE_TYPE.MCLC]: ({ a, b, c, alpha, beta, gamma }) => {
        const cosAlpha = math.cos((alpha / 180) * math.PI);
        return [
            [a / 2, b / 2, 0],
            [-a / 2, b / 2, 0],
            [0, c * cosAlpha, c * math.sqrt(1 - cosAlpha * cosAlpha)],
        ];
    },

    // Algorithm from http://pymatgen.org/_modules/pymatgen/core/lattice.html (from_params)
    [LATTICE_TYPE.TRI]: ({ a, b, c, alpha, beta, gamma }) => {
        // convert angles to Radians
        [alpha, beta, gamma] = [alpha, beta, gamma].map(
            (x) => math.unit(x, "degree").to("rad").value,
        );

        const [cosAlpha, cosBeta, cosGamma] = [alpha, beta, gamma].map((x) => math.cos(x));
        const [sinAlpha, sinBeta, sinGamma] = [alpha, beta, gamma].map((x) => math.sin(x));

        const gammaStar = math.acos((cosAlpha * cosBeta - cosGamma) / (sinAlpha * sinBeta));
        const cosGammaStar = math.cos(gammaStar);
        const sinGammaStar = math.sin(gammaStar);

        return [
            [a * sinBeta, 0.0, a * cosBeta],
            [-b * sinAlpha * cosGammaStar, b * sinAlpha * sinGammaStar, b * cosAlpha],
            [0.0, 0.0, c],
        ];
    },

    // alternative implementation
    [`${LATTICE_TYPE.TRI}alt`]: ({ a, b, c, alpha, beta, gamma }) => {
        const cosAlpha = math.cos((alpha / 180) * math.PI);
        const cosBeta = math.cos((beta / 180) * math.PI);
        const cosGamma = math.cos((gamma / 180) * math.PI);
        const sinGamma = math.sqrt(1 - cosGamma * cosGamma);
        return [
            [a, 0.0, 0.0],
            [b * cosGamma, b * sinGamma, 0.0],
            [
                c * cosBeta,
                (c / sinGamma) * (cosAlpha - cosBeta * cosGamma),
                (c / sinGamma) *
                    math.sqrt(
                        sinGamma * sinGamma -
                            cosAlpha * cosAlpha -
                            cosBeta * cosBeta +
                            2 * cosAlpha * cosBeta * cosGamma,
                    ),
            ],
        ];
    },
};

/**
 * Returns lattice vectors for a primitive cell for a lattice.
 * @param {Lattice} lattice - Lattice instance.
 * @param {Boolean[]} skipRounding - whether to skip rounding the lattice vectors.
 * @return {Array[]} Cell.vectorsAsArray
 */
export function primitiveCell(lattice, skipRounding = false) {
    const [vectorA, vectorB, vectorC] = PRIMITIVE_CELLS[lattice.type || LATTICE_TYPE.TRI](lattice);
    // set precision and remove JS floating point artifacts
    !skipRounding &&
        [vectorA, vectorB, vectorC].map((vec) =>
            vec.map((c) => math.precise(c)).map(math.roundToZero),
        );
    return [vectorA, vectorB, vectorC];
}
