import { LatticeImplicitSchema, LatticeTypeSchema } from "@mat3ra/esse/dist/js/types";

import constants from "../constants";
import math from "../math";
import { LATTICE_TYPE_CONFIGS, Vector, VectorsAsArray } from "./types";

export type Units = Required<LatticeImplicitSchema>["units"];

export interface FromVectorsProps {
    a: Vector; // vector of the lattice.
    b: Vector; // vector of the lattice.
    c: Vector; // vector of the lattice.
    alat?: number; // scaling factor for the vector coordinates, defaults to 1.
    units?: Units["length"]; // units for the coordinates (eg. angstrom, crystal).
    type?: LatticeTypeSchema; // type of the lattice to be created, defaults to TRI when not provided.
    skipRounding?: boolean; // whether to skip rounding the resulting lattice values, defaults to `false`.
}

/*
 * @summary: class that holds parameters of a Bravais Lattice: a, b, c, alpha, beta, gamma + corresponding units.
 * When stored as class variables units for lengths are always "angstrom"s, angle - "degree"s
 */
export class LatticeBravais implements LatticeImplicitSchema {
    a: number;

    b: number;

    c: number;

    alpha: number;

    beta: number;

    gamma: number;

    units: Units;

    type: LatticeImplicitSchema["type"];

    /**
     * Create a Bravais lattice.
     */
    constructor(config: Partial<LatticeImplicitSchema>) {
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

        this.a = a * k;
        this.b = b * k;
        this.c = c * k;
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;
        this.type = type;
        this.units = units;
    }

    static _roundValue(x: number): number {
        return math.precise(math.roundToZero(x));
    }

    /**
     * Create a Bravais lattice from vectors.
     */
    static fromVectors({
        a,
        b,
        c,
        alat = 1,
        units = "angstrom",
        type = "TRI",
        skipRounding = false,
    }: FromVectorsProps) {
        const roundValue = skipRounding ? (x: number) => x : this._roundValue;
        const config = {
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
                angle: "degree" as Units["angle"],
            },
        };
        return new (this.prototype.constructor as typeof LatticeBravais)(config);
    }

    /**
     * See fromVectors above.
     */
    static fromVectorArrays(array: VectorsAsArray, type: LatticeTypeSchema, skipRounding = true) {
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
        )?.editables;
        // ["a", "gamma"] => {a: true, gamma: true}
        if (editablesList) {
            editablesList.forEach((element) => {
                Object.assign(object, {
                    [element]: true,
                });
            });
        }

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
    toJSON(): LatticeImplicitSchema {
        return {
            ...this,
            units: {
                length: this.units.length,
                angle: this.units.angle,
            },
        };
    }
}
