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

const DEFAULT_CELL: [PointSchema, PointSchema, PointSchema] = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
];

export class Cell implements CellSchema {
    alat = 1;

    private _a: RoundedVector3D;

    private _b: RoundedVector3D;

    private _c: RoundedVector3D;

    constructor(
        a: PointSchema = DEFAULT_CELL[0],
        b: PointSchema = DEFAULT_CELL[1],
        c: PointSchema = DEFAULT_CELL[2],
    ) {
        this._a = new RoundedVector3D(a);
        this._b = new RoundedVector3D(b);
        this._c = new RoundedVector3D(c);
    }

    static fromVectorsArray(vectors: PointSchema[]): Cell {
        return new Cell(vectors[0], vectors[1], vectors[2]);
    }

    get a(): PointSchema {
        return this._a.value_rounded;
    }

    get b(): PointSchema {
        return this._b.value_rounded;
    }

    get c(): PointSchema {
        return this._c.value_rounded;
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
        return math.roundArrayOrNumber(this.volume, 9) as number;
    }

    clone(): Cell {
        return Cell.fromVectorsArray(this.vectorArrays);
    }

    cloneAndScaleByMatrix(matrix: number[][]) {
        const newCell = this.clone();
        newCell.scaleByMatrix(matrix);
        return newCell;
    }

    /**
     * Convert a point (in crystal coordinates) to cartesian.
     */
    convertPointToCartesian(point: PointSchema): CodeMath.MathType {
        return MULT(point, this.vectorArrays);
    }

    /**
     * Convert a point (in cartesian coordinates) to crystal (fractional).
     */
    convertPointToFractional(point: PointSchema): PointSchema {
        return MULT(point, INV(this.vectorArrays)) as unknown as PointSchema;
    }

    /**
     * Check whether a point is inside the cell.
     * @param point - the point to conduct the check for.
     * @param tolerance - numerical tolerance.
     */
    isPointInsideCell(point: PointSchema, tolerance = constants.tolerance.length): boolean {
        return (
            this.convertPointToFractional(point)
                .map((c: number) => math.isBetweenZeroInclusiveAndOne(c, tolerance))
                // @ts-ignore
                .reduce((a: boolean, b: boolean): boolean => a && b)
        );
    }

    /**
     * Returns the index of the cell vector, most collinear with the testVector.
     * @param testVector
     */
    getMostCollinearVectorIndex(testVector: Vector): number {
        const angles = this.vectorArrays.map((v) => math.angleUpTo90(v, testVector, "deg"));
        return angles.findIndex((el: number) => el === math.min(angles));
    }

    /**
     * Scale this cell by right-multiplying it to a matrix (nested array)
     */
    scaleByMatrix(matrix: number[][]) {
        // @ts-ignore
        [this.a, this.b, this.c] = MATRIX_MULT(matrix, this.vectorArrays);
    }
}
