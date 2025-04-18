import { RoundedVector3D } from "@mat3ra/code";
import {
    LatticeExplicitUnit,
    LatticeImplicitSchema,
    LatticeTypeEnum,
} from "@mat3ra/esse/dist/js/types";

import { Cell } from "../cell/cell";
import { primitiveCell } from "../cell/primitive_cell";
import constants from "../constants";
import math from "../math";
import { LatticeVectorUnits, Vector, VectorsAsArray } from "../types";

export interface LatticeBravaisConfig extends Required<LatticeImplicitSchema> {
    isConventional?: boolean;
}

export type LatticeVectorsConfig = LatticeExplicitUnit & {
    type?: LatticeTypeEnum;
};

export class LatticeVectors extends Cell {
    units: LatticeVectorUnits;

    constructor(config: LatticeExplicitUnit) {
        super(config.a, config.b, config.c);
        const { a, b, c, alat = 1, units = "angstrom" } = config;
        const k = constants.units.bohr === units ? constants.coefficients.BOHR_TO_ANGSTROM : 1;

        this._a = new RoundedVector3D(a.map((x) => x * k));
        this._b = new RoundedVector3D(b.map((x) => x * k));
        this._c = new RoundedVector3D(c.map((x) => x * k));
        this.alat = alat;
        this.units = "angstrom"; // always store as angstrom internally
    }

    static _roundValue(arr: number[]): number[] {
        return arr.map((el) => math.precise(math.roundToZero(el)));
    }

    /**
     * Construct a LatticeVectors object from Bravais parameters.
     * Optionally generate primitive or conventional cell.
     */
    static fromBravais({
        a = 1,
        b = a,
        c = a,
        alpha = 90,
        beta = alpha,
        gamma = alpha,
        units = {
            length: "angstrom",
            angle: "degree",
        },
        type = "TRI",
        isConventional = false,
    }: LatticeBravaisConfig): LatticeVectors {
        // fall back to primitive if conventional not explicitly requested
        if (!isConventional) type = "TRI";

        const [vecA, vecB, vecC] = primitiveCell({
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
            a: vecA,
            b: vecB,
            c: vecC,
            alat: 1,
        });
    }

    /**
     * Return raw vector values as array of 3-vectors
     */
    get vectorArrays(): VectorsAsArray {
        return [this.a, this.b, this.c];
    }

    /**
     * Serialize with rounded coordinates
     */
    toJSON(): Required<LatticeExplicitUnit> {
        return {
            a: LatticeVectors._roundValue(this.a) as Vector,
            b: LatticeVectors._roundValue(this.b) as Vector,
            c: LatticeVectors._roundValue(this.c) as Vector,
            alat: this.alat,
            units: this.units,
        };
    }
}
