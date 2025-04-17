import {
    LatticeExplicitUnit,
    LatticeImplicitSchema,
    LatticeTypeEnum,
} from "@mat3ra/esse/dist/js/types";

import { primitiveCell } from "../cell/primitive_cell";
import constants from "../constants";
import math from "../math";
import { latticeVectorUnits, Vector, VectorsAsArray } from "../types";

export interface LatticeBravaisConfig extends Required<LatticeImplicitSchema> {
    isConventional?: boolean;
}

export type LatticeVectorsConfig = LatticeExplicitUnit & {
    type?: LatticeTypeEnum;
};

export class LatticeVectors implements Required<LatticeExplicitUnit> {
    a: Vector;

    b: Vector;

    c: Vector;

    alat: number;

    units: latticeVectorUnits;

    constructor(config: LatticeExplicitUnit) {
        const { a, b, c, alat = 1, units = "angstrom" } = config;
        const k = constants.units.bohr === units ? constants.coefficients.BOHR_TO_ANGSTROM : 1;

        this.a = a.map((x) => x * k) as Vector;
        this.b = b.map((x) => x * k) as Vector;
        this.c = c.map((x) => x * k) as Vector;
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
