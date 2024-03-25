import { LatticeExplicitUnit, LatticeImplicitSchema } from "@mat3ra/esse/dist/js/types";

import { primitiveCell } from "../cell/primitive_cell";
import constants from "../constants";
import math from "../math";
import { Vector } from "./types";

type RequiredLatticeExplicitUnit = Required<LatticeExplicitUnit>;

export interface BravaisConfigProps extends Partial<LatticeImplicitSchema> {
    isConventional?: boolean;
}

/*
 * @summary: class that holds parameters of a Bravais Lattice: a, b, c, alpha, beta, gamma + corresponding units.
 * When stored as class variables units for lengths are always "angstrom"s, angle - "degree"s
 */
export class LatticeVectors implements RequiredLatticeExplicitUnit {
    a: RequiredLatticeExplicitUnit["a"];

    b: RequiredLatticeExplicitUnit["b"];

    c: RequiredLatticeExplicitUnit["c"];

    alat: RequiredLatticeExplicitUnit["alat"];

    units: RequiredLatticeExplicitUnit["units"];

    /**
     * Create a Bravais lattice.
     */
    constructor(config: LatticeExplicitUnit) {
        const { a, b, c, alat = 1, units = "angstrom" } = config;
        const k = constants.units.bohr === units ? constants.coefficients.BOHR_TO_ANGSTROM : 1;

        this.a = a.map((x) => x * k) as Vector;
        this.b = b.map((x) => x * k) as Vector;
        this.c = c.map((x) => x * k) as Vector;
        this.alat = alat;
        this.units = "angstrom";
    }

    static _roundValue(arr: number[]): number[] {
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
        type = "TRI",
        isConventional = false,
    }: BravaisConfigProps) {
        // use "direct" lattice constructor for primitive lattice
        // eslint-disable-next-line no-param-reassign
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

        return new LatticeVectors({
            a: vectorA,
            b: vectorB,
            c: vectorC,
            alat: 1,
        });
    }

    get vectorArrays(): [Vector, Vector, Vector] {
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
    toJSON(): RequiredLatticeExplicitUnit {
        return {
            ...this,
            a: LatticeVectors._roundValue(this.a),
            b: LatticeVectors._roundValue(this.b),
            c: LatticeVectors._roundValue(this.c),
        };
    }
}
