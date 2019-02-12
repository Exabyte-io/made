import math from "../math";

const ADD = math.add;

/**
 * Counts integer shifts in both positive and negative directions by sweeping the [-10, 10] interval and
 * obtaining the boundaries for the sub-interval where shifted point located on the lattice represented by
 * latticeVectors is within the cell.
 * TODO: implement an optimized version and auto-locate the amplitude.
 * @param {Cell} cell - initial cell.
 * @param {Array} point - point in 3D.
 * @param {Cell} anotherCell - supercell.
 * @return {Array[]} - Nested array of integer shifts.
 * @example [[1, 0, 0], [-1, 0, 0]]
 */
function generateTranslationCombinations(cell, point, anotherCell, amplitude = 10) {

    const combinations = [];
    const range = Array.from({length: 2 * amplitude + 1}, (v, k) => k - amplitude);

    range.forEach(i => {
        range.forEach(j => {
            range.forEach(k => {
                const shift = cell.convertPointToCartesian([i, j, k]);
                if (anotherCell.isPointInsideCell(ADD(point, shift))) {
                    combinations.push([i, j, k]);
                }
            })
        })
    });

    return combinations;
}

export default {
    generateTranslationCombinations,
}
