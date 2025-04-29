import { RoundedVector3D } from "@mat3ra/code";
import { math } from "@mat3ra/code/dist/js/math";
import {
    Coordinate3DSchema,
    LatticeVectorsSchema as CellSchema,
    Matrix3X3Schema,
    Vector3DSchema,
} from "@mat3ra/esse/dist/js/types";

import constants from "../constants";

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

    static fromVectorsArray(vectors: Matrix3X3Schema): Cell {
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

    get vectorArrays(): Matrix3X3Schema {
        return [this._a.value, this._b.value, this._c.value];
    }

    get vectorArraysRounded(): Matrix3X3Schema {
        return [this._a.valueRounded, this._b.valueRounded, this._c.valueRounded];
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

    cloneAndScaleByMatrix(matrix: Matrix3X3Schema) {
        const newCell = this.clone();
        newCell.scaleByMatrix(matrix);
        return newCell;
    }

    convertPointToCartesian(point: Coordinate3DSchema): Coordinate3DSchema {
        return MULT(point, this.vectorArrays) as unknown as Coordinate3DSchema;
    }

    convertPointToCrystal(point: Coordinate3DSchema): Coordinate3DSchema {
        return MULT(point, INV(this.vectorArrays)) as unknown as Coordinate3DSchema;
    }

    isPointInsideCellCartesian(point: Coordinate3DSchema): boolean {
        const { tolerance } = this.constructor as typeof Cell;
        return (
            this.convertPointToCrystal(point)
                .map((c: number) => math.isBetweenZeroInclusiveAndOne(c, tolerance))
                // @ts-ignore
                .reduce((a: boolean, b: boolean): boolean => a && b)
        );
    }

    // eslint-disable-next-line class-methods-use-this
    isPointInsideCellCrystal(point: Coordinate3DSchema): boolean {
        const { tolerance } = this.constructor as typeof Cell;
        return (
            point
                .map((c: number) => math.isBetweenZeroInclusiveAndOne(c, tolerance))
                // @ts-ignore
                .reduce((a: boolean, b: boolean): boolean => a && b)
        );
    }

    isPointInsideCell(point: Coordinate3DSchema, useCrystal = false): boolean {
        if (useCrystal) {
            return this.isPointInsideCellCrystal(point);
        }
        return this.isPointInsideCellCartesian(point);
    }

    getMostCollinearVectorIndex(testVector: Vector3DSchema): number {
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
