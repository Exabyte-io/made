import CodeMath, { math } from "@mat3ra/code/dist/js/math";

import { Coordinate } from "../basis/types";
import constants from "../constants";
import { Vector, VectorsAsArray } from "../lattice/types";

const MATRIX = math.matrix;
const MULT = math.multiply;
const INV = math.inv;
// @ts-ignore
const MATRIX_MULT = (...args) => MULT(...args.map((x) => MATRIX(x))).toArray();

type Point = Coordinate | CodeMath.Matrix | CodeMath.MathType;

/*
 * Cell represents a unit cell in geometrical form: 3x3 matrix where rows are cell vectors.
 * Example: [[1, 0, 0], [0, 1, 0], [0, 0, 1]].
 */
export class Cell {
    tolerance = 1;

    vector1: Vector;

    vector2: Vector;

    vector3: Vector;

    /**
     * Create a cell.
     * @param nestedArray {Number[][]} is an array of cell vectors in cartesian Angstrom units.
     */
    constructor(nestedArray: VectorsAsArray) {
        [this.vector1, this.vector2, this.vector3] = nestedArray;
    }

    /**
     * Get cell vectors as (a nested) array.
     * @example [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
     */
    get vectorsAsArray(): VectorsAsArray {
        return [this.vector1, this.vector2, this.vector3] as VectorsAsArray;
    }

    clone(): Cell {
        return new (this.constructor as typeof Cell)(this.vectorsAsArray);
    }

    cloneAndScaleByMatrix(matrix: number[][]) {
        const newCell = this.clone();
        newCell.scaleByMatrix(matrix);
        return newCell;
    }

    /**
     * Convert a point (in crystal coordinates) to cartesian.
     */
    convertPointToCartesian(point: Point): CodeMath.MathType {
        return MULT(point, this.vectorsAsArray);
    }

    /**
     * Convert a point (in cartesian coordinates) to crystal (fractional).
     */
    convertPointToFractional(point: Point): Coordinate {
        return MULT(point, INV(this.vectorsAsArray)) as Coordinate;
    }

    /**
     * Check whether a point is inside the cell.
     * @param point - the point to conduct the check for.
     * @param tolerance - numerical tolerance.
     */
    isPointInsideCell(point: Point, tolerance = constants.tolerance.length): boolean {
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
        // @ts-ignore
        const angles = this.vectorsAsArray.map((v) => math.angleUpTo90(v, testVector));
        return angles.findIndex((el) => el === math.min(angles));
    }

    /**
     * Scale this cell by right-multiplying it to a matrix (nested array)
     */
    scaleByMatrix(matrix: number[][]) {
        // @ts-ignore
        [this.vector1, this.vector2, this.vector3] = MATRIX_MULT(matrix, this.vectorsAsArray);
    }
}
