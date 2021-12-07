import constants from "../constants";
import math from "../math";
import { LATTICE_TYPE_CONFIGS } from "./types";

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
        const editablesList = LATTICE_TYPE_CONFIGS.find((entry) => entry.code === this.type)
            .editables;
        // ["a", "gamma"] => {a: true, gamma: true}
        editablesList.forEach((element) => {
            Object.assign(object, {
                [element]: true,
            });
        });
        return object;
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
