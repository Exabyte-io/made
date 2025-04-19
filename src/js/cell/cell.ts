import { RoundedVector3D } from "@mat3ra/code";
import CodeMath, { math } from "@mat3ra/code/dist/js/math";
import { LatticeExplicitUnit as CellSchema, PointSchema } from "@mat3ra/esse/dist/js/types";

import constants from "../constants";
import { Vector } from "../types";

const MATRIX = math.matrix;
const MULT = math.multiply;
const INV = math.inv;
// @ts-ignore
const MATRIX_MULT = (...args) => MULT(...args.map((x) => MATRIX(x))).toArray();

export class Cell implements CellSchema {
    static RoundedVector3DClassReference = RoundedVector3D;

    static roundPrecision = 9;

    static tolerance = constants.tolerance.length; // in crystal coordinates

    static defaultConfig: CellSchema = {
        a: [1.0, 0, 0],
        b: [0, 1.0, 0],
        c: [0, 0, 1.0],
    };

    a: CellSchema["a"];

    b: CellSchema["b"];

    c: CellSchema["c"];

    alat = 1;

    units: CellSchema["units"] = "angstrom";

    constructor(config = Cell.defaultConfig) {
        const { a, b, c } = config;
        this.a = a;
        this.b = b;
        this.c = c;
    }

    static fromVectorsArray(vectors: PointSchema[]): Cell {
        return new Cell({
            a: vectors[0],
            b: vectors[1],
            c: vectors[2],
        });
    }

    get _a() {
        return new Cell.RoundedVector3DClassReference(this.a);
    }

    get _b() {
        return new Cell.RoundedVector3DClassReference(this.b);
    }

    get _c() {
        return new Cell.RoundedVector3DClassReference(this.c);
    }

    get vectorArrays(): PointSchema[] {
        return [this._a.value, this._b.value, this._c.value];
    }

    get vectorArraysRounded(): PointSchema[] {
        return [this._a.value_rounded, this._b.value_rounded, this._c.value_rounded];
    }

    get volume(): number {
        return math.det(this.vectorArrays);
    }

    get volumeRounded(): number {
        return math.roundArrayOrNumber(this.volume, Cell.roundPrecision) as number;
    }

    clone(): Cell {
        return new Cell({ a: this.a, b: this.b, c: this.c });
    }

    cloneAndScaleByMatrix(matrix: number[][]) {
        const newCell = this.clone();
        newCell.scaleByMatrix(matrix);
        return newCell;
    }

    convertPointToCartesian(point: PointSchema): CodeMath.MathType {
        return MULT(point, this.vectorArrays);
    }

    convertPointToCrystal(point: PointSchema): PointSchema {
        return MULT(point, INV(this.vectorArrays)) as unknown as PointSchema;
    }

    isPointInsideCellCartesian(point: PointSchema): boolean {
        const { tolerance } = this.constructor as typeof Cell;
        return (
            this.convertPointToCrystal(point)
                .map((c: number) => math.isBetweenZeroInclusiveAndOne(c, tolerance))
                // @ts-ignore
                .reduce((a: boolean, b: boolean): boolean => a && b)
        );
    }

    // eslint-disable-next-line class-methods-use-this
    isPointInsideCellCrystal(point: PointSchema): boolean {
        const { tolerance } = this.constructor as typeof Cell;
        return (
            point
                .map((c: number) => math.isBetweenZeroInclusiveAndOne(c, tolerance))
                // @ts-ignore
                .reduce((a: boolean, b: boolean): boolean => a && b)
        );
    }

    isPointInsideCell(point: PointSchema, useCrystal = false): boolean {
        if (useCrystal) {
            return this.isPointInsideCellCrystal(point);
        }
        return this.isPointInsideCellCartesian(point);
    }

    getMostCollinearVectorIndex(testVector: Vector): number {
        const angles = this.vectorArrays.map((v) => math.angleUpTo90(v, testVector, "deg"));
        return angles.findIndex((el: number) => el === math.min(angles));
    }

    scaleByMatrix(matrix: number[][]) {
        [this.a, this.b, this.c] = MATRIX_MULT(matrix, this.vectorArrays);
    }

    toJSON(): CellSchema {
        return {
            a: this._a.toJSON(),
            b: this._b.toJSON(),
            c: this._c.toJSON(),
            alat: this.alat,
            units: this.units,
        };
    }
}
