import constants from "../constants";
import math from "../math";
import { LATTICE_TYPE, LATTICE_TYPE_CONFIGS } from "./types";

/*
 * @summary: class that holds parameters of a Bravais Lattice: a, b, c, alpha, beta, gamma + corresponding units.
 * When stored as class variables units for lengths are always "angstrom"s, angle - "degree"s
 */
export class LatticeBravais {
    /**
     * Create a Bravais lattice.
     * @param {Object} config - Config object.
     * @param {Number} config.a - lattice constant a.
     * @param {Number} config.b - lattice constant b.
     * @param {Number} config.c - lattice constant c.
     * @param {Number} config.alpha - lattice angle alpha.
     * @param {Number} config.beta - lattice angle beta.
     * @param {Number} config.gamma - lattice angle gamma.
     * @param {Object} config.units - units container.
     * @param {String} config.units.length - units for length (eg. "angstrom", "crystal").
     * @param {String} config.units.angle - units for angles (eg. "degree", "radian").
     */
    constructor(config) {
        const {
            a = 1, // default lattice is cubic with unity in edge sizes
            b = a,
            c = a,
            alpha = 90,
            beta = alpha,
            gamma = alpha,
            // if we do not know what lattice type this is => set to TRI
            type = "TRI",
            units = {
                length: "angstrom",
                angle: "degree",
            },
        } = config;
        const k =
            constants.units.bohr === units.length ? constants.coefficients.BOHR_TO_ANGSTROM : 1;
        Object.assign(this, {
            a: a * k,
            b: b * k,
            c: c * k,
            alpha,
            beta,
            gamma,
            type,
            units,
        });
    }

    static _roundValue(x) {
        return math.precise(math.roundToZero(x));
    }

    /**
     * Create a Bravais lattice from vectors.
     * @param {Object} - Config object.
     * @param {Array} a - vector of the lattice.
     * @param {Array} b - vector of the lattice.
     * @param {Array} c - vector of the lattice.
     * @param {Number} alat - scaling factor for the vector coordinates, defaults to 1.
     * @param {String} units - units for the coordinates (eg. angstrom, crystal).
     * @param {String} type - type of the lattice to be created, defaults to TRI when not provided.
     * @param {Boolean} skipRounding - whether to skip rounding the resulting lattice values, defaults to `false`.
     */
    static fromVectors({
        a,
        b,
        c,
        alat = 1,
        units = "angstrom",
        type = "TRI",
        skipRounding = false,
    }) {
        const roundValue = skipRounding ? (x) => x : this._roundValue;
        return new this.prototype.constructor({
            a: roundValue(math.vlen(a) * alat),
            b: roundValue(math.vlen(b) * alat),
            c: roundValue(math.vlen(c) * alat),
            alpha: roundValue(math.angle(b, c, "deg")),
            beta: roundValue(math.angle(a, c, "deg")),
            gamma: roundValue(math.angle(a, b, "deg")),
            // initially we do not know what lattice type this is => set to TRI
            type,
            units: {
                length: units,
                angle: "degree",
            },
        });
    }

    /**
     * See fromVectors above.
     */
    static fromVectorArrays(array, type, skipRounding = true) {
        return this.fromVectors({
            a: array[0],
            b: array[1],
            c: array[2],
            type,
            skipRounding, // do not round the values to avoid loosing precision by default
        });
    }

    /**
     * Get the list of editable keys (eg. 'a', 'alpha') for the current lattice.
     * @return {Object}
     * @example {a: true, b: false, c: false, alpha: true, beta: false, gamma: false}
     */
    get editables() {
        const object = {};
        const editablesList = LATTICE_TYPE_CONFIGS.find(
            (entry) => entry.code === this.type,
        ).editables;
        // ["a", "gamma"] => {a: true, gamma: true}
        editablesList.forEach((element) => {
            Object.assign(object, {
                [element]: true,
            });
        });
        return object;
    }

    /**
     * Return lattice type from vectors array.
     * @param {Number[][]} array
     * @returns {String}
     */
    // eslint-disable-next-line no-unused-vars
    static typeFromArrays(array) {
        const type = "unknown";
        // TODO: implement
        return type;
    }

    /**
     * Return primitive lattice vectors from type and necessary parameters of: a, b, c, cosbc, cosac, cosab.
     * According to this: https://arxiv.org/pdf/1004.2974.pdf
     * @param {string} type
     * @param {number} a
     * @param {number} [b]
     * @param {number} [c]
     * @param {number} [alpha] - angle between b and c in radians
     * @param {number} [beta] - angle between a and c in radians
     * @param {number} [gamma] - angle between a and b in radians
     * @returns {{vectors: number[][], type: string}}
     */
    // eslint-disable-next-line no-unused-vars
    static vectorsFromType(type, a, b, c, alpha, beta, gamma) {
        let vectors = [];
        // TODO: check if this many math.js operations give correct result

        switch (type) {
            case LATTICE_TYPE.CUB:
                vectors = [
                    [a, 0, 0],
                    [0, a, 0],
                    [0, 0, a],
                ];
                break;
            case LATTICE_TYPE.FCC:
                vectors = [
                    [0, a / 2, a / 2],
                    [a / 2, 0, a / 2],
                    [a / 2, a / 2, 0],
                ];
                break;
            case LATTICE_TYPE.TET:
                vectors = [
                    [a, 0, 0],
                    [0, a, 0],
                    [0, 0, c],
                ];
                break;
            case LATTICE_TYPE.BCT:
                vectors = [
                    [-a / 2, a / 2, c / 2],
                    [a / 2, -a / 2, c / 2],
                    [a / 2, a / 2, -c / 2],
                ];
                break;
            case LATTICE_TYPE.ORC:
                vectors = [
                    [a, 0, 0],
                    [0, b, 0],
                    [0, 0, c],
                ];
                break;
            case LATTICE_TYPE.ORCF:
                vectors = [
                    [0, b / 2, c / 2],
                    [a / 2, 0, c / 2],
                    [a / 2, b / 2, 0],
                ];
                break;
            case LATTICE_TYPE.ORCI:
                vectors = [
                    [-a / 2, b / 2, c / 2],
                    [a / 2, -b / 2, c / 2],
                    [a / 2, b / 2, -c / 2],
                ];
                break;
            case LATTICE_TYPE.ORCC:
                vectors = [
                    [a / 2, -b / 2, 0],
                    [a / 2, b / 2, 0],
                    [0, 0, c],
                ];
                break;
            case LATTICE_TYPE.HEX:
                vectors = [
                    [a / 2, (-a * math.sqrt(3)) / 2, 0],
                    [a / 2, (a * math.sqrt(3)) / 2, 0],
                    [0, 0, c],
                ];
                break;
            case LATTICE_TYPE.RHL:
                vectors = [
                    [a * math.cos(alpha / 2), -a * math.sin(alpha / 2), 0],
                    [a * math.cos(alpha / 2), a * math.sin(alpha / 2), 0],
                    [
                        (a * math.cos(alpha)) / math.cos(alpha / 2),
                        0,
                        a * math.sqrt(1 - math.cos(alpha) ** 2 / math.cos(alpha / 2) ** 2),
                    ],
                ];
                break;
            case LATTICE_TYPE.MCL:
                vectors = [
                    [a, 0, 0],
                    [0, b, 0],
                    [0, c * math.cos(alpha), c * math.sin(alpha)],
                ];
                break;
            case LATTICE_TYPE.MCLC:
                vectors = [
                    [a / 2, b / 2, 0],
                    [-a / 2, b / 2, 0],
                    [0, c * math.cos(alpha), c * math.sin(alpha)],
                ];

                break;
            case LATTICE_TYPE.TRI:
                // eslint-disable-next-line no-case-declarations
                const z =
                    (c / math.sin(gamma)) *
                    math.sqrt(
                        math.sin(gamma) ** 2 -
                            math.cos(alpha) ** 2 -
                            math.cos(beta) ** 2 +
                            2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma),
                    );
                vectors = [
                    [a, 0, 0],
                    [b * math.cos(gamma), b * math.sin(gamma), 0],
                    [
                        c * math.cos(beta),
                        (c / math.sin(gamma)) *
                            (math.cos(alpha) - math.cos(beta) * math.cos(gamma)),
                        z,
                    ],
                ];
                break;
            default:
                throw new Error("Unknown lattice type");
        }
        return { vectors, type };
    }

    /**
     * Serialize class instance to JSON.
     * @example As below:
         {
            "a" : 3.867,
            "b" : 3.867,
            "c" : 3.867,
            "alpha" : 60,
            "beta" : 60,
            "gamma" : 60,
            "units" : {
                "length" : "angstrom",
                "angle" : "degree"
            },
            "type" : "FCC"
         }
     */
    toJSON() {
        return {
            ...this,
            units: {
                length: this.units.length,
                angle: this.units.angle,
            },
        };
    }
}
