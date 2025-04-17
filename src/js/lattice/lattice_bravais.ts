import {
    LatticeImplicitSchema as LatticeBravaisSchema,
    LatticeTypeEnum,
} from "@mat3ra/esse/dist/js/types";

import constants from "../constants";
import math from "../math";
import { LatticeVectorsConfig, VectorsAsArray } from "./lattice_vectors";
import { LATTICE_TYPE_CONFIGS } from "./lattice_types";

export type Units = Required<LatticeBravaisSchema>["units"];

/*
 * @summary: class that holds parameters of a Bravais Lattice: a, b, c, alpha, beta, gamma + corresponding units.
 * When stored as class variables units for lengths are always "angstrom"s, angle - "degree"s
 */
export class LatticeBravais implements LatticeBravaisSchema {
    a: number;

    b: number;

    c: number;

    alpha: number;

    beta: number;

    gamma: number;

    units: Units;

    type: LatticeBravaisSchema["type"];

    /**
     * Create a Bravais lattice.
     */
    constructor(config: LatticeBravaisSchema) {
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
    }: LatticeVectorsConfig): LatticeBravais {
        const round = this._roundValue;

        const config: LatticeBravaisSchema = {
            a: round(math.vlen(a) * alat),
            b: round(math.vlen(b) * alat),
            c: round(math.vlen(c) * alat),
            alpha: round(math.angle(b, c, "deg")),
            beta: round(math.angle(a, c, "deg")),
            gamma: round(math.angle(a, b, "deg")),
            type: type ?? "TRI",
            units: {
                length: units,
                angle: "degree",
            },
        };

        return new LatticeBravais(config);
    }

    /**
     * See fromVectors above.
     */
    static fromVectorArrays(array: VectorsAsArray, type: LatticeTypeEnum) {
        return this.fromVectors({
            a: array[0],
            b: array[1],
            c: array[2],
            type,
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
            (entry: any) => entry.code === this.type,
        )?.editables;
        // ["a", "gamma"] => {a: true, gamma: true}
        if (editablesList) {
            editablesList.forEach((element: any) => {
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
    toJSON(): LatticeBravaisSchema {
        return {
            ...this,
            units: {
                length: this.units.length,
                angle: this.units.angle,
            },
        };
    }
}
