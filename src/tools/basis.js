import _ from "underscore";
import {Basis} from "../basis/basis";
import math from "../math";

const ADD = math.add;
const MULT = math.multiply;

/**
 * Returns a repeated basis of a crystal.
 * @param basis {Basis} Original basis.
 * @param repetitions{Number[]} Repetition vector `[x, y, z]`, in each spatial dimension.
 * @return {Basis} New Basis.
 */
function repeat(basis, repetitions) {
    let i, j, k, shiftI = 0, shiftJ = 0, shiftK = 0;

    // clone original basis and assert it is in cartesian coordinates
    const newBasis = basis.clone();
    const basisCloneInCrystalCoordinates = basis.clone();

    newBasis.toCrystal();
    basisCloneInCrystalCoordinates.toCrystal();

    for (i = 1; i <= repetitions[0]; i++) {
        for (j = 1; j <= repetitions[1]; j++) {
            for (k = 1; k <= repetitions[2]; k++) {
                // for each atom in original basis add one with a repetition
                basisCloneInCrystalCoordinates.elements.forEach((element, index) => {
                    const coord = basisCloneInCrystalCoordinates.getCoordinateByIndex(index);
                    // only add atoms if shifts are non-zero
                    (shiftI || shiftJ || shiftK) && newBasis.addAtom({
                        element: element,
                        coordinate: [
                            coord[0] + shiftI,
                            coord[1] + shiftJ,
                            coord[2] + shiftK
                        ]
                    });
                });
                shiftK += 1;
            }
            shiftK = 0;
            shiftJ += 1;
        }
        shiftJ = 0;
        shiftI += 1;
    }

    if (basis.isInCartesianUnits) newBasis.toCartesian();

    return newBasis;
}

/**
 * Returns a set of Bases for a crystal interpolated from initial to final crystal.
 *          Can be used to generate atomic configurations along a chemical reaction path, for example.
 * @param initialBasis {Basis} Original initialBasis.
 * @param finalBasis {Basis} Final initialBasis.
 * @param numberOfSteps{Number} Number of intermediate steps.
 * @return {Basis[]} List of all bases.
 */
function interpolate(initialBasis, finalBasis, numberOfSteps = 1) {
    // check that initial and final basis have the same cell
    if (!initialBasis.hasEquivalentCellTo(finalBasis))
        throw new Error("basis.interpolate: Basis cells are not equal");

    // clone original initialBasis and assert it is in cartesian coordinates
    const initialBasisCopy = initialBasis.clone();
    const finalBasisCopy = finalBasis.clone();

    initialBasisCopy.toCrystal();
    finalBasisCopy.toCrystal();

    const initialCoordinates = _.flatten(initialBasisCopy.coordinatesAsArray);
    const finalCoordinates = _.flatten(finalBasisCopy.coordinatesAsArray);
    const delta = ADD(finalCoordinates, MULT(initialCoordinates, -1));

    const resultingListOfBases = [];

    for (let i = 1; i <= numberOfSteps; i++) {
        const normalizedStepIndex = i / (numberOfSteps + 1);
        const intermediateCoordinates = _linearInterpolation(initialCoordinates, delta, normalizedStepIndex);
        const vectorSize = 3;
        const intermediateCoordinatesAsNestedArray = _.toArray(
            _.groupBy(intermediateCoordinates, (a, b) => Math.floor(b / vectorSize))
        );
        const intermediateBasis = initialBasis.clone();
        intermediateBasis.coordinates = intermediateCoordinatesAsNestedArray;
        resultingListOfBases.push(intermediateBasis)
    }

    return resultingListOfBases;
}

/**
 * Calculates linear function `y = kx + b` for vectors. Isolated for modularity.
 * @param initialCoordinates {Array} - b.
 * @param delta {Array} - x.
 * @param normalizedStepIndex {Number} - k.
 * @return {Basis[]} List of all bases.
 */

function _linearInterpolation(initialCoordinates, delta, normalizedStepIndex) {
    return ADD(initialCoordinates, MULT(delta, normalizedStepIndex))
}

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
 * @param basisCoordinates
 * @returns {*[]}
 */
export function centerOfCoordinates(basisCoordinates) {
    const transposedBasisCoordinates = math.transpose(basisCoordinates);
    const centerOfCoordinatesVectors = [];
    for (let i = 0; i < transposedBasisCoordinates.length; i++) {
        let center = transposedBasisCoordinates[i].reduce((a, b) => a + b) / basisCoordinates.length;
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
 * @param basisCoordinates
 */
export function pairwiseDistance(basisCoordinates) {
    const maxDistance = 0;
    for (let i = 0; i < basisCoordinates.length; i++) {
        for (let j = i + 1; j < basisCoordinates.length; j++) {
            const distance = math.vDist(basisCoordinates[i], basisCoordinates[j]);
            if (distance > maxDistance) {
                let maxDistance = distance;
            }
        }
    }
    return maxDistance;
}

export default {
    repeat,
    interpolate,
    centerOfCoordinates,
    pairwiseDistance
}
