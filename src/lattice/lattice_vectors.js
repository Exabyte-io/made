import constants from "../constants";
import {primitiveCell} from "../cell/primitive_cell";
import math from "../math";

/*
 * @summary: class that holds parameters of a Bravais Lattice: a, b, c, alpha, beta, gamma + corresponding units.
 * When stored as class variables units for lengths are always "angstrom"s, angle - "degree"s
 */
export class LatticeVectors {
    constructor(config) {
        const {a, b, c, alat = 1, units = "angstrom"} = config;
        const k = constants.units.bohr === units ? constants.coefficients.BOHR_TO_ANGSTROM : 1;
        Object.assign(this, {
            a: a.map(c => c * k),
            b: b.map(c => c * k),
            c: c.map(c => c * k),
            alat,
            units: "angstrom",
        });
    }

    static _roundValue(arr) {
        return arr.map(el => math.precise(math.roundToZero(el)));
    }

    /*
     * @summary Constructs a Bravais lattice from lattice vectors
     *          Supports conventional lattice as parameters and primitive too.
     *          Algorithm from http://pymatgen.org/_modules/pymatgen/core/lattice.html (from_params)
     */
    static fromBravais({
                           a = 1,  // default lattice is cubic with unity in edge sizes
                           b = a,
                           c = a,
                           alpha = 90,
                           beta = 90,
                           gamma = 90,
                           units = {
                               length: "angstrom",
                               angle: "degree"
                           },
                           type,
                           isConventional = false,
                       }) {

        // use "direct" lattice constructor for primitive lattice
        if (!isConventional) type = 'TRI';

        // set precision and remove JS floating point artifacts
        const [vectorA, vectorB, vectorC] = primitiveCell({
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
            units,
            type
        });

        return new LatticeVectors.prototype.constructor({
            a: vectorA,
            b: vectorB,
            c: vectorC,
            alat: 1,
            unis: units.length,
        });
    }

    get vectorArrays() {
        return [
            this.a,
            this.b,
            this.c,
        ]
    }

    toJSON() {
        return Object.assign({}, this, {
            a: LatticeVectors._roundValue(this.a),
            b: LatticeVectors._roundValue(this.b),
            c: LatticeVectors._roundValue(this.c),
        });
    }
}
