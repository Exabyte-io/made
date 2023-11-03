import constants from "../constants";
import math from "../math";
import { LATTICE_TYPE_CONFIGS, LatticeType, Vector, VectorsAsArray } from "./types";

export interface BravaisConfig {
    a?: number; // lattice constant a.
    b?: number; // lattice constant b.
    c?: number; // lattice constant c.
    alpha?: number; // lattice angle alpha.
    beta?: number; // lattice angle beta.
    gamma?: number; // lattice angle gamma.
    units?: {
        // units container
        length: string; // units for length (eg. "angstrom", "crystal")
        angle: string; // units for angles (eg. "degree", "radian")
    };
    type?: LatticeType;
}

export type RequiredBravaisConfig = Required<BravaisConfig>;

export interface FromVectorsProps {
    a: Vector; // vector of the lattice.
    b: Vector; // vector of the lattice.
    c: Vector; // vector of the lattice.
    alat?: number; // scaling factor for the vector coordinates, defaults to 1.
    units?: string; // units for the coordinates (eg. angstrom, crystal).
    type?: LatticeType; // type of the lattice to be created, defaults to TRI when not provided.
    skipRounding?: boolean; // whether to skip rounding the resulting lattice values, defaults to `false`.
}

/*
 * @summary: class that holds parameters of a Bravais Lattice: a, b, c, alpha, beta, gamma + corresponding units.
 * When stored as class variables units for lengths are always "angstrom"s, angle - "degree"s
 */
export class LatticeBravais implements RequiredBravaisConfig {
    a: number;

    b: number;

    c: number;

    alpha: number;

    beta: number;

    gamma: number;

    units: {
        length: string;
        angle: string;
    };

    type: LatticeType;

    /**
     * Create a Bravais lattice.
     */
    constructor(config: BravaisConfig) {
        const {
            a = 1, // default lattice is cubic with unity in edge sizes
            b = a,
            c = a,
            alpha = 90,
            beta = alpha,
            gamma = alpha,
            // if we do not know what lattice type this is => set to TRI
            type = LatticeType.TRI,
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
        type = LatticeType.TRI,
        skipRounding = false,
    }: FromVectorsProps) {
        const roundValue = skipRounding ? (x: number) => x : this._roundValue;

        return new (this.prototype.constructor as typeof LatticeBravais)({
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
    static fromVectorArrays(array: VectorsAsArray, type: LatticeType, skipRounding = true) {
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
    toJSON(): RequiredBravaisConfig {
        return {
            ...this,
            units: {
                length: this.units.length,
                angle: this.units.angle,
            },
        };
    }
}
