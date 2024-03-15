"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
const underscore_1 = __importDefault(require("underscore"));
// eslint-disable-next-line no-unused-vars
const basis_1 = require("../basis/basis");
const math_1 = __importDefault(require("../math"));
const ADD = math_1.default.add;
const MULT = math_1.default.multiply;
/**
 * Returns a repeated basis of a crystal.
 * @param basis {Basis} Original basis.
 * @param repetitions{Number[]} Repetition vector `[x, y, z]`, in each spatial dimension.
 * @return {Basis} New Basis.
 */
function repeat(basis, repetitions) {
    let i, j, k;
    let shiftI = 0;
    let shiftJ = 0;
    let shiftK = 0;
    // clone original basis and assert it is in cartesian coordinates
    const newBasis = basis.clone();
    const basisCloneInCrystalCoordinates = basis.clone();
    newBasis.toCrystal();
    basisCloneInCrystalCoordinates.toCrystal();
    for (i = 1; i <= repetitions[0]; i += 1) {
        for (j = 1; j <= repetitions[1]; j += 1) {
            for (k = 1; k <= repetitions[2]; k += 1) {
                // for each atom in original basis add one with a repetition
                // eslint-disable-next-line no-loop-func
                basisCloneInCrystalCoordinates.elements.forEach((element, index) => {
                    const coord = basisCloneInCrystalCoordinates.getCoordinateByIndex(index);
                    // only add atoms if shifts are non-zero
                    if (shiftI || shiftJ || shiftK) {
                        newBasis.addAtom({
                            element,
                            coordinate: [coord[0] + shiftI, coord[1] + shiftJ, coord[2] + shiftK],
                        });
                    }
                });
                shiftK += 1;
            }
            shiftK = 0;
            shiftJ += 1;
        }
        shiftJ = 0;
        shiftI += 1;
    }
    if (basis.isInCartesianUnits)
        newBasis.toCartesian();
    return newBasis;
}
/**
 * Calculates linear function `y = kx + b` for vectors. Isolated for modularity.
 * @param initialCoordinates {Array} - b.
 * @param delta {Array} - x.
 * @param normalizedStepIndex {Number} - k.
 * @return {Basis[]} List of all bases.
 */
function _linearInterpolation(initialCoordinates, delta, normalizedStepIndex) {
    return ADD(initialCoordinates, MULT(delta, normalizedStepIndex));
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
    const initialCoordinates = underscore_1.default.flatten(initialBasisCopy.coordinatesAsArray);
    const finalCoordinates = underscore_1.default.flatten(finalBasisCopy.coordinatesAsArray);
    const delta = ADD(finalCoordinates, MULT(initialCoordinates, -1));
    const resultingListOfBases = [];
    for (let i = 1; i <= numberOfSteps; i++) {
        const normalizedStepIndex = i / (numberOfSteps + 1);
        const intermediateCoordinates = _linearInterpolation(initialCoordinates, delta, normalizedStepIndex);
        const vectorSize = 3;
        const intermediateCoordinatesAsNestedArray = underscore_1.default.toArray(underscore_1.default.groupBy(intermediateCoordinates, (a, b) => Math.floor(b / vectorSize)));
        const intermediateBasis = initialBasis.clone();
        intermediateBasis.coordinates = intermediateCoordinatesAsNestedArray;
        resultingListOfBases.push(intermediateBasis);
    }
    return resultingListOfBases;
}
exports.default = {
    repeat,
    interpolate,
};
//# sourceMappingURL=basis.js.map