import constants from "../constants";
import math from "../math";
import {LATTICE_TYPE_CONFIGS} from "./types";

/*
 * @summary: class that holds parameters of a Bravais Lattice: a, b, c, alpha, beta, gamma + corresponding units.
 * When stored as class variables units for lengths are always "angstrom"s, angle - "degree"s
 */
export class LatticeBravais {
    constructor(config) {
        const {
            a = 1, // default lattice is cubic with unity in edge sizes
            b = a,
            c = a,
            alpha = 90,
            beta = alpha,
            gamma = alpha,
            // if we do not know what lattice type this is => set to TRI
            type = 'TRI',
            units = {
                length: 'angstrom',
                angle: 'degree'
            }
        } = config;
        const k = constants.units.bohr === units.length ? constants.coefficients.BOHR_TO_ANGSTROM : 1;
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

    static _roundValue (x) {
        return math.precise(math.roundToZero(x));
    }

    static fromVectors({a, b, c, alat = 1, units = 'angstrom', type = 'TRI', skipRounding = false}) {
        const roundValue = skipRounding ? x => x : this._roundValue;
        return new this.prototype.constructor({
            a: roundValue(math.vlen(a) * alat),
            b: roundValue(math.vlen(b) * alat),
            c: roundValue(math.vlen(c) * alat),
            alpha: roundValue(math.angle(b, c, 'deg')),
            beta: roundValue(math.angle(a, c, 'deg')),
            gamma: roundValue(math.angle(a, b, 'deg')),
            // initially we do not know what lattice type this is => set to TRI
            type,
            units: {
                length: units,
                angle: 'degree'
            }
        });
    }

    static fromVectorArrays(array, type, skipRounding = true) {
        return this.fromVectors({
            a: array[0],
            b: array[1],
            c: array[2],
            type,
            skipRounding,  // do not round the values to avoid loosing precision by default
        });
    }

    get editables() {
        const object = {};
        const editablesList = LATTICE_TYPE_CONFIGS.find(entry => entry.code === this.type).editables;
        // ["a", "gamma"] => {a: true, gamma: true}
        editablesList.forEach(element => {
            Object.assign(object, {
                [element]: true,
            })
        });
        return object;
    }

    toJSON() {
        return Object.assign({}, this, {
            units: {
                length: this.units.length,
                angle: this.units.angle,
            }
        })
    }
}
