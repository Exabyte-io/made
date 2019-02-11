import math from "../math";
import constants from "../constants";

const MATRIX = math.matrix;
const MULT = math.multiply;
const INV = math.inv;
const ADD = math.add;
const MATRIX_MULT = (...args) => MULT(...args.map(x => MATRIX(x))).toArray();

/*
 * @summary Cell represents a unit cell in geometrical form: 3x3 matrix where rows are cell vectors:
 *          Example: [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
 *
 * @param nestedArray {Number[][]} is an array of cell vectors in cartesian Angstrom units
 */
export class Cell {

    tolerance = 1;

    constructor(nestedArray) {
        this.vector1 = nestedArray[0];
        this.vector2 = nestedArray[1];
        this.vector3 = nestedArray[2];
    }

    get vectorsAsArray() {
        return [
            this.vector1,
            this.vector2,
            this.vector3
            // assert that no near-zero artifacts are present (ie. 1.6 x 10^-16) before attempting inversion
            // @param tolerance {Number} Maximum tolerance to small numbers, used to avoid artifacts on matrix inversion
        ].map(v => v.map(c => (math.abs(c) < constants.tolerance.lengthAngstrom) ? 0 : c));
    }

    clone() {return new this.constructor(this.vectorsAsArray)}

    cloneAndScaleByMatrix(matrix) {
        const newCell = this.clone();
        newCell.scaleByMatrix(matrix);
        return newCell;
    }

    convertPointToCartesian(point) {
        return MULT(point, this.vectorsAsArray);
    }

    convertPointToFractional(point) {
        return MULT(point, INV(this.vectorsAsArray));
    }

    isPointInsideCell(point, tolerance = constants.tolerance.length) {
        return this.convertPointToFractional(point).map(c => math.isBetweenZeroInclusiveAndOne(c, tolerance))
            .reduce((a, b) => a && b);
    }

    /*
     * @summary Counts integer shifts in both positive and negative directions by sweeping the [-10, 10] interval and
     *          obtaining the boundaries for the sub-interval where shifted point located on the lattice represented by
     *          latticeVectors is within the cell.
     *          TODO: implement an optimized version and auto-locate amplitude.
     *          TODO: migrate to tools
     */
    generateTranslationCombinations(point, anotherCell, amplitude = 10) {

        const combinations = [];
        const range = Array.from({length: 2 * amplitude + 1}, (v, k) => k - amplitude);

        range.forEach(i => {
            range.forEach(j => {
                range.forEach(k => {
                    const shift = this.convertPointToCartesian([i, j, k]);
                    if (anotherCell.isPointInsideCell(ADD(point, shift))) {
                        combinations.push([i, j, k]);
                    }
                })
            })
        });

        return combinations;
    }

    getMostCollinearVectorIndex(testVector) {
        const angles = this.vectorsAsArray.map(v => math.angleUpTo90(v, testVector));
        return angles.findIndex(el => el === math.min(angles));
    }

    /*
     * @summary Scaled this cell by right-multiplying it to a matrix (nested array)
     */
    scaleByMatrix(matrix) {
        [this.vector1, this.vector2, this.vector3] = MATRIX_MULT(matrix, this.vectorsAsArray);
    }

}
