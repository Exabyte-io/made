import { primitiveCell } from "../cell/primitive_cell";
import constants from "../constants";
import math from "../math";
import { BravaisConfig } from "./lattice_bravais";
import { LatticeType, Vector } from "./types";

interface LatticeVectorsConfig {
    a: Vector; // vector of the lattice.
    b: Vector; // vector of the lattice.
    c: Vector; // vector of the lattice.
    alat: number; // scaling factor for the vector coordinates, defaults to 1.
    units?: string; // units container.
}

export type RequiredLatticeVectorsConfig = Required<LatticeVectorsConfig>;

export interface BravaisConfigProps extends BravaisConfig {
    isConventional?: boolean;
}

/*
 * @summary: class that holds parameters of a Bravais Lattice: a, b, c, alpha, beta, gamma + corresponding units.
 * When stored as class variables units for lengths are always "angstrom"s, angle - "degree"s
 */
export class LatticeVectors implements LatticeVectorsConfig {
    a: Vector;

    b: Vector;

    c: Vector;

    alat: number;

    units: string;

    /**
     * Create a Bravais lattice.
     */
    constructor(config: LatticeVectorsConfig) {
        const { a, b, c, alat = 1, units = "angstrom" } = config;
        const k = constants.units.bohr === units ? constants.coefficients.BOHR_TO_ANGSTROM : 1;

        this.a = a.map((x) => x * k) as Vector;
        this.b = b.map((x) => x * k) as Vector;
        this.c = c.map((x) => x * k) as Vector;
        this.alat = alat;
        this.units = "angstrom";
    }

    static _roundValue(arr: number[]) {
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
        type = LatticeType.TRI,
        isConventional = false,
    }: BravaisConfigProps) {
        // use "direct" lattice constructor for primitive lattice
        // eslint-disable-next-line no-param-reassign
        if (!isConventional) type = LatticeType.TRI;

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
    toJSON(): RequiredLatticeVectorsConfig {
        return {
            ...this,
            a: LatticeVectors._roundValue(this.a),
            b: LatticeVectors._roundValue(this.b),
            c: LatticeVectors._roundValue(this.c),
        };
    }
}
