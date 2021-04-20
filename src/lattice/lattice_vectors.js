import { primitiveCell } from "../cell/primitive_cell";
import constants from "../constants";
import math from "../math";

/*
 * @summary: class that holds parameters of a Bravais Lattice: a, b, c, alpha, beta, gamma + corresponding units.
 * When stored as class variables units for lengths are always "angstrom"s, angle - "degree"s
 */
export class LatticeVectors {
    /**
     * Create a Bravais lattice.
     * @param {Object} config - Config object.
     * @param {Array} config.a - vector of the lattice.
     * @param {Array} config.b - vector of the lattice.
     * @param {Array} config.c - vector of the lattice.
     * @param {Number} config.alat - scaling factor for the vector coordinates, defaults to 1.
     * @param {String} config.units - units container.
     */
    constructor(config) {
        const { a, b, c, alat = 1, units = "angstrom" } = config;
        const k = constants.units.bohr === units ? constants.coefficients.BOHR_TO_ANGSTROM : 1;
        Object.assign(this, {
            a: a.map((x) => x * k),
            b: b.map((x) => x * k),
            c: c.map((x) => x * k),
            alat,
            units: "angstrom",
        });
    }

    static _roundValue(arr) {
        return arr.map((el) => math.precise(math.roundToZero(el)));
    }

    /*
     * Constructs a Bravais lattice from lattice vectors
     * Supports conventional lattice as parameters and primitive too.
     * Algorithm from http://pymatgen.org/_modules/pymatgen/core/lattice.html (from_params)
     * For parameters see `LatticeBravais.constructor`.
     */
    static fromBravais({
        a = 1, // default lattice is cubic with unity in edge sizes
        b = a,
        c = a,
        alpha = 90,
        beta = 90,
        gamma = 90,
        units = {
            length: "angstrom",
            angle: "degree",
        },
        type,
        isConventional = false,
    }) {
        // use "direct" lattice constructor for primitive lattice
        if (!isConventional) type = "TRI";

        // set precision and remove JS floating point artifacts
        const [vectorA, vectorB, vectorC] = primitiveCell({
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
            units,
            type,
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
        return [this.a, this.b, this.c];
    }

    /**
     * Serialize class instance to JSON.
     * @example As below:
     {
            "a" : [
                3.34892,
                0,
                1.9335
            ],
            "b" : [
                1.116307,
                3.157392,
                1.9335
            ],
            "c" : [
                0,
                0,
                3.867
            ],
            "alat" : 1,
            "units" : "angstrom"
        }
     */
    toJSON() {
        return {
            ...this,
            a: LatticeVectors._roundValue(this.a),
            b: LatticeVectors._roundValue(this.b),
            c: LatticeVectors._roundValue(this.c),
        };
    }
}
