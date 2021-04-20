import { tolerance as tol } from "../constants";
import math from "../math";

const MATRIX = math.matrix;
const MULT = math.multiply;
const INV = math.inv;
const MATRIX_MULT = (...args) => MULT(...args.map((x) => MATRIX(x))).toArray();

/*
 * Cell represents a unit cell in geometrical form: 3x3 matrix where rows are cell vectors.
 * Example: [[1, 0, 0], [0, 1, 0], [0, 0, 1]].
 */
export class Cell {
    tolerance = 1;

    /**
     * Create a cell.
     * @param nestedArray {Number[][]} is an array of cell vectors in cartesian Angstrom units.
     */
    constructor(nestedArray) {
        [this.vector1, this.vector2, this.vector3] = nestedArray;
    }

    /**
     * Get cell vectors as (a nested) array.
     * @return {Array[]}
     * @example [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
     */
    get vectorsAsArray() {
        return [
            this.vector1,
            this.vector2,
            this.vector3,
            // assert that no near-zero artifacts are present (ie. 1.6 x 10^-16) before attempting inversion
            // @param tolerance {Number} Maximum tolerance to small numbers, used to avoid artifacts on matrix inversion
        ].map((v) => v.map((c) => (math.abs(c) < tol.lengthAngstrom ? 0 : c)));
    }

    clone() {
        return new this.constructor(this.vectorsAsArray);
    }

    cloneAndScaleByMatrix(matrix) {
        const newCell = this.clone();
        newCell.scaleByMatrix(matrix);
        return newCell;
    }

    /**
     * Convert a point (in crystal coordinates) to cartesian.
     * @return {Array}
     */
    convertPointToCartesian(point) {
        return MULT(point, this.vectorsAsArray);
    }

    /**
     * Convert a point (in cartesian coordinates) to crystal (fractional).
     * @return {Array}
     */
    convertPointToFractional(point) {
        return MULT(point, INV(this.vectorsAsArray));
    }

    /**
     * Check whether a point is inside the cell.
     * @param {Array} point - the point to conduct the check for.
     * @param {Number} tolerance - numerical tolerance.
     * @return {Boolean}
     */
    isPointInsideCell(point, tolerance = tol.length) {
        return this.convertPointToFractional(point)
            .map((c) => math.isBetweenZeroInclusiveAndOne(c, tolerance))
            .reduce((a, b) => a && b);
    }

    /**
     * Returns the index of the cell vector, most collinear with the testVector.
     * @param testVector
     * @return {Number}
     */
    getMostCollinearVectorIndex(testVector) {
        const angles = this.vectorsAsArray.map((v) => math.angleUpTo90(v, testVector));
        return angles.findIndex((el) => el === math.min(angles));
    }

    /**
     * Scale this cell by right-multiplying it to a matrix (nested array)
     */
    scaleByMatrix(matrix) {
        [this.vector1, this.vector2, this.vector3] = MATRIX_MULT(matrix, this.vectorsAsArray);
    }
}
