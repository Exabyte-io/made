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
 * @summary Function takes initial basis coordinates and transposes them so that all X, Y, and Z dimension
 * values are in their own nested array. Then the sums of each nested array are computed and used to calculate
 * the center point shift needed to move a molecule so that it is centered around a point (0,0,0)
 *
 * initialCoordinates
 * [[x1, y1, z1],
 *  [x2, y2, z2],
 *  [.., .., ..],
 *  [xn, yn, zn]]
 *
 * transposedInitialCoordaintes
 * [[x1, x2, ...xn],
 *  [y1, y2, ...yn],
 *  [z1, z2, ...zn]]
 *
 * center point x = sum(all x elements) / no. x elements
 *
 * @param basis
 * @returns {*}
 */
function centerPointShift(basis) {
    const transposedInitialCoordinates = math.transpose(basis);
    const centerPointShift = [];
    for (let i = 0; i < transposedInitialCoordinates.length; i++) {
        let center = transposedInitialCoordinates[i].reduce((a, b) => a + b) / basis.length;
        centerPointShift.push(center);
    }
    return centerPointShift
}

/**
 * @summary function returns the max distance between pairs of basis elements.
 * @param basis
 */
function maxPairwiseDistance(basis) {
    const maxDistance = 0;
    for (let i = 0; i < len(this._coordinates); i++) {
        for (let j = 0; j < len(this._coordinates); j++) {
            const distance = math.vDist(this._coordinates[i], this._coordinates[j]);
            if (distance >  maxDistance) {
                let maxDistance = distance;
            }
        }
    }
    return maxDistance
}

export default {
    repeat,
    interpolate,
    centerPointShift,
    maxPairwiseDistance
}
