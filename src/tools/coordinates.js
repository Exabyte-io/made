import math from "../math";

/**
 * @summary Function takes basis coordinates and transposes them so that the values for each dimension of the
 *  the basis are in their own nested array. Then the center point for each dimension of the coordinates is calculated.
 *
 * initial basisCoordinates
 * [[x1, y1, z1],
 *  [x2, y2, z2],
 *  [.., .., ..],
 *  [xn, yn, zn]]
 *
 * transposed basisCoordinates
 * [[x1, x2, ...xn],
 *  [y1, y2, ...yn],
 *  [z1, z2, ...zn]]
 *
 * Returns an array = [xCenter, yCenter, zCenter]
 *
 * @param {Array} basisCoordinatesArray
 * @param {Number} nAtoms
 * @returns {Array}
 */
export function coordinatesGetCenterOfSpaceAsVector(basisCoordinatesArray, nAtoms) {
    const coordinates = []
    const transposedBasisCoordinates = math.transpose(basisCoordinatesArray);
    const centerOfCoordinatesVectors = [];
    for (let i = 0; i < 3; i++) {
        let center = transposedBasisCoordinates[i].reduce((a, b) => a + b) / nAtoms;
        centerOfCoordinatesVectors.push(center);
    }
    return centerOfCoordinatesVectors;
}

/**
 * @summary function returns the max distance between pairs of basis coordinates
 * basis coordinates = [[x1, y1, z1], [x2, y2, z2], ... [xn, yn, zn]]
 *
 * calculates distance between each basis coordinate set.
 * n = last set of coordinates
 * n-1 = second to last set of coordinates
 *
 * Iterate through pairs without redundancy.
 *      pair 0,1   pair 0,2  pair 0,3 ... pair 0,n
 *      -          pair 1,2  pair 1,3 ... pair 1,n
 *      -     -    -         pair 2,3 ... pair 2,n
 *      -     -    -         ...      ...
 *      -     -    -         ...      ... pair n-1, n
 *
 * @param {ArrayWithIds} basisCoordinates
 * @param {Number} nAtoms
 * @return {Number}
 */
export function coordinatesGetMaxPairwiseDistance(basisCoordinates, nAtoms) {
    let maxDistance = 0;
    for (let i = 0; i < nAtoms; i++) {
        for (let j = i + 1; j < nAtoms; j++) {
            const distance = math.vDist(basisCoordinates.getArrayElementByIndex(i), basisCoordinates.getArrayElementByIndex(j));
            if (distance > maxDistance) {
                maxDistance = distance;
            }
        }
    }
    return maxDistance;
}

export default {
    coordinatesGetCenterOfSpaceAsVector,
    coordinatesGetMaxPairwiseDistance
}
