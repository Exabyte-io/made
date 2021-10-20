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
 * @summary finds the displacement needed for the molecule to be centered around (0,0,0)
 * @param initialCoordinates
 * @returns {number[]}
 */
function centeredMolecule(initialCoordinates) {
    const nAtoms = initialCoordinates.length;
    const mass = 1.0;
    let atomCounter = 0;
    let xSum = 0;
    let ySum = 0;
    let zSum = 0;
    while (atomCounter < nAtoms) {
        xSum += initialCoordinates[atomCounter[0]]
        ySum += initialCoordinates[atomCounter[1]]
        zSum += initialCoordinates[atomCounter[2]]
        atomCounter++;
    }

    const xCenter = x_sum / (mass * nAtoms);
    const yCenter = y_sum / (mass * nAtoms);
    const zCenter = z_sum / (mass * nAtoms);
    const centerOfStructureShift = [xCenter, yCenter, zCenter];

    return centerOfStructureShift
}

function moleculeMaxRadius(initialCoordinates) {
    let maxRadius = 0.;
    let atomCounterA = 0;
    let atomCounterB = 0;
    const nAtoms = initialCoordinates.length;
    while (atomCounterA < nAtoms) {
        while (atomCounterB < nAtoms) {
            const atomA = initialCoordinates[atomCounterA];
            const atomB = initialCoordinates[atomCounterB];
            const distance = math.distance(atomA, atomB);
            if (distance > maxRadius) {
                maxRadius = distance;
            }
            atomCounterB++;
        }
        atomCounterA++;
    }
    maxRadius = math.ceil(maxRadius);
    return maxRadius
}

export default {
    repeat,
    interpolate,
    moleculeMaxRadius,
    centeredMolecule,
}
