"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.Cell = void 0;
const math_1 = require("@mat3ra/code/dist/js/math");
const constants_1 = __importDefault(require("../constants"));
const MATRIX = math_1.math.matrix;
const MULT = math_1.math.multiply;
const INV = math_1.math.inv;
// @ts-ignore
const MATRIX_MULT = (...args) => MULT(...args.map((x) => MATRIX(x))).toArray();
/*
 * Cell represents a unit cell in geometrical form: 3x3 matrix where rows are cell vectors.
 * Example: [[1, 0, 0], [0, 1, 0], [0, 0, 1]].
 */
class Cell {
    /**
     * Create a cell.
     * @param nestedArray {Number[][]} is an array of cell vectors in cartesian Angstrom units.
     */
    constructor(nestedArray) {
        this.tolerance = 1;
        [this.vector1, this.vector2, this.vector3] = nestedArray;
    }
    /**
     * Get cell vectors as (a nested) array.
     * @example [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
     */
    get vectorsAsArray() {
        return [
            this.vector1,
            this.vector2,
            this.vector3,
            // assert that no near-zero artifacts are present (ie. 1.6 x 10^-16) before attempting inversion
            // @param tolerance {Number} Maximum tolerance to small numbers, used to avoid artifacts on matrix inversion
        ].map((v) => v.map((c) => (math_1.math.abs(c) < constants_1.default.tolerance.lengthAngstrom ? 0 : c)));
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
     */
    convertPointToCartesian(point) {
        return MULT(point, this.vectorsAsArray);
    }
    /**
     * Convert a point (in cartesian coordinates) to crystal (fractional).
     */
    convertPointToFractional(point) {
        return MULT(point, INV(this.vectorsAsArray));
    }
    /**
     * Check whether a point is inside the cell.
     * @param point - the point to conduct the check for.
     * @param tolerance - numerical tolerance.
     */
    isPointInsideCell(point, tolerance = constants_1.default.tolerance.length) {
        return (this.convertPointToFractional(point)
            .map((c) => math_1.math.isBetweenZeroInclusiveAndOne(c, tolerance))
            // @ts-ignore
            .reduce((a, b) => a && b));
    }
    /**
     * Returns the index of the cell vector, most collinear with the testVector.
     * @param testVector
     */
    getMostCollinearVectorIndex(testVector) {
        // @ts-ignore
        const angles = this.vectorsAsArray.map((v) => math_1.math.angleUpTo90(v, testVector));
        return angles.findIndex((el) => el === math_1.math.min(angles));
    }
    /**
     * Scale this cell by right-multiplying it to a matrix (nested array)
     */
    scaleByMatrix(matrix) {
        // @ts-ignore
        [this.vector1, this.vector2, this.vector3] = MATRIX_MULT(matrix, this.vectorsAsArray);
    }
}
exports.Cell = Cell;
