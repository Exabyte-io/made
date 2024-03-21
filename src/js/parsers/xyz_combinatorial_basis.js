/* eslint-disable max-classes-per-file */
import _ from "underscore";
import * as s from "underscore.string";

import { Basis } from "../basis/basis";
import math from "../math";

/**
 * @summary Combinatorial XYZ basis class and related. Create and get all information about basis and elements in it.
 * Constructor accepts string in extended XYZ format. Extended XYZ format is as follows:
 *
 * 1. Regular XYZ
 * ```
 *      Si 0.0 0.0 0.0
 *      Li 0.5 0.5 0.5
 * ```
 * 2. Permutation (slash-separated)
 * ```
 *      Si/Ge/As 0.0 0.0 0.0
 *      Si/Ge    0.5 0.5 0.0
 * ```
 * 3. Combination (comma-separated)
 * ```
 *      Si,Ge,As 0.0 0.0 0.0
 *      Si,Ge    0.5 0.5 0.0
 * ```
 * More at: https://exabyte.atlassian.net/wiki/display/PD/Combinatorial+materials+set
 */

/**
 * Regular expression for basis line.
 * @type {RegExp}
 */
// eslint-disable-next-line max-len
const LINE_REGEX =
    /^([A-Z][a-z]?\/?,?)+\s+(-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)\s*$/gi;
// vacancy characters will be used to create vacancies on basis generation
const VACANCY_CHARACTER = "VAC";
const COMBINATION_DELIMITER = ",";
const PERMUTATION_DELIMITER = "/";
/**
 * Basis validation error codes.
 * @type {{MIXING_IN_SINGLE_LINE: number, MIXING_IN_MULTI_LINES: number, REGEX_NOT_PASSED: number}}
 */
const ERROR_CODES = {
    MIXING_IN_SINGLE_LINE: 1,
    MIXING_IN_MULTI_LINES: 2,
    REGEX_NOT_PASSED: 3,
};

export class WrongBasisFormat extends Error {
    constructor(xyz, message, type) {
        super(message, type);
        this.xyz = xyz;
    }
}

export class CombinatorialBasis {
    /**
     * Creates Combinatorial basis
     * @param eXYZ
     */
    constructor(eXYZ) {
        this._xyz = eXYZ;
        this._lines = s
            .lines(eXYZ)
            .map((x) => x.trim())
            .filter((x) => x !== "")
            .map(this._parseBasisLine);

        this._hasPermutationLine = this._lines.reduce((mem, a) => {
            return mem || a.isPermutation;
        }, false);
        this._hasCombinationLine = this._lines.reduce((mem, a) => {
            return mem || a.isCombination;
        }, false);

        if (this._hasPermutationLine && this._hasCombinationLine) {
            throw new WrongBasisFormat(
                this._xyz,
                "Basis contains mixed permutation and combination.",
                ERROR_CODES.MIXING_IN_MULTI_LINES,
            );
        }
    }

    /**
     * Parses combinatorial basis line and returns result as Object. Throws exception if line is not valid:
     * - does not meet RegExp
     * - mixes permutation and combination
     * @param str {String} Combinatorial basis' line
     * @param index {Number} order of the
     * @return {{displayName: string, isCombination: boolean, isPermutation: boolean, elements: Array, coordinates: *[]}}
     * @private
     */
    _parseBasisLine(str, index) {
        if (!str.match(LINE_REGEX)) {
            throw new WrongBasisFormat(
                this._xyz,
                `Line #${index + 1}: "${str}" contains errors. ` +
                    'Allowed formats: "Si 0 0 0", "Si/Li 0.5 0.5 0.5", "Si,Ge 0.7 0.7 0.8"',
                ERROR_CODES.REGEX_NOT_PASSED,
            );
        }

        const containsPermutation = str.indexOf(PERMUTATION_DELIMITER) > -1;
        const containsCombination = str.indexOf(COMBINATION_DELIMITER) > -1;

        if (containsCombination && containsPermutation) {
            throw new WrongBasisFormat(
                this._xyz,
                `Line #${index} contains mixed permutation and combination.`,
                ERROR_CODES.MIXING_IN_SINGLE_LINE,
            );
        }

        let elements = [];

        const words = s
            .words(str)
            .map((x) => x.trim())
            .filter((x) => x != null);

        if (containsCombination) {
            elements = words[0].split(COMBINATION_DELIMITER).map((x) => x.trim());
        } else if (containsPermutation) {
            elements = words[0].split(PERMUTATION_DELIMITER).map((x) => x.trim());
        } else {
            elements = [words[0]];
        }

        const coordinates = [parseFloat(words[1]), parseFloat(words[2]), parseFloat(words[3])];

        return {
            displayName: `ELEMENT_${index}`,
            isCombination: containsCombination,
            isPermutation: containsPermutation,
            elements,
            coordinates,
        };
    }

    /**
     * Returns array of ALL unique elements used in basis.
     * @return {String[]}
     */
    get uniqueElements() {
        return _.chain(this._lines)
            .map((line) => line.elements)
            .flatten()
            .unique()
            .value()
            .sort();
    }

    static toBasisConfig(array, units = "crystal", cell = Basis.defaultCell) {
        return {
            elements: _.pluck(array, "element"),
            coordinates: _.pluck(array, "coordinates"),
            units,
            cell,
        };
    }

    /**
     * Returns array of regular bases extracted from current combinatorial basis.
     * @return {Basis[]|Object[]}
     */
    get allBasisConfigs() {
        let result = [];
        if (this._hasPermutationLine) {
            result = this._permutation();
        } else if (this._hasCombinationLine) {
            result = this._combination();
        } else {
            const items = [];
            this._lines.forEach((line) => {
                items.push({
                    element: line.elements[0],
                    coordinates: line.coordinates,
                });
            });
            result = [items];
        }
        return result.map((x) => CombinatorialBasis.toBasisConfig(x));
    }

    /**
     * Returns array of regular bases extracted from current combinatorial basis with combinations.
     * @private
     */
    _combination() {
        const dimensions = [];
        this._lines.forEach((line) => {
            const itemsSet = [];
            line.elements.forEach((element) => {
                // omit vacancy characters
                itemsSet.push({
                    element,
                    coordinates: line.coordinates,
                });
            });
            dimensions.push(itemsSet);
        });
        const basisSet = math.cartesianProduct.apply(null, dimensions);
        return basisSet.map((basis) =>
            basis.filter((entry) => entry.element !== VACANCY_CHARACTER),
        );
    }

    /**
     * Returns array of regular bases extracted from current combinatorial basis with permutations.
     * @return {Basis[]}
     * @private
     */
    _permutation() {
        const maxLen = Math.max(...this._lines.map((x) => x.elements.length));
        const bases = [];
        for (let i = 0; i < maxLen; i++) {
            const items = [];
            this._lines.forEach((line) => {
                const element =
                    line.elements.length <= i ? _.last(line.elements) : line.elements[i];
                if (element !== VACANCY_CHARACTER) {
                    items.push({
                        element,
                        coordinates: line.coordinates,
                    });
                }
            });
            bases.push(items);
        }
        return bases;
    }

    /**
     * Returns true if current combinatorial basis contains more than one regular basis.
     * @return {Boolean}
     */
    isCombinatorial() {
        return this._hasCombinationLine || this._hasPermutationLine;
    }
}
