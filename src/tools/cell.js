import math from "../math";
import {tolerance as TOLERANCE} from "../constants";

/**
 * @summary Returns true if cells are equal within tolerance, otherwise - false.
 * @param cellVectors1 {Array} Cell Vectors 1, eg. [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
 * @param cellVectors2 {Array} Cell Vectors 2.
 * @param tolerance{Number} Tolerance
 * @return {Boolean}
 */
function areEqual(cellVectors1, cellVectors2, tolerance = TOLERANCE.pointsDistance) {
    const equalWithTolerance = (vec1, vec2) => (math.vDist(vec1, vec2) <= tolerance);
    return !cellVectors1.map((vector, idx) => (equalWithTolerance(vector, cellVectors2[idx])))
        .some(x => !x);
}

export default {
    areEqual,
}
