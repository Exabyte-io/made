export class WrongBasisFormat extends Error {
    constructor(xyz: any, message: any, type: any);
    xyz: any;
}
export class CombinatorialBasis {
    static toBasisConfig(array: any, units?: string, cell?: [import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema]): {
        elements: any[];
        coordinates: any[];
        units: string;
        cell: [import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema];
    };
    /**
     * Creates Combinatorial basis
     * @param eXYZ
     */
    constructor(eXYZ: any);
    _xyz: any;
    _lines: {
        displayName: string;
        isCombination: boolean;
        isPermutation: boolean;
        elements: any[];
        coordinates: any[];
    }[];
    _hasPermutationLine: boolean;
    _hasCombinationLine: boolean;
    /**
     * Parses combinatorial basis line and returns result as Object. Throws exception if line is not valid:
     * - does not meet RegExp
     * - mixes permutation and combination
     * @param str {String} Combinatorial basis' line
     * @param index {Number} order of the
     * @return {{displayName: string, isCombination: boolean, isPermutation: boolean, elements: Array, coordinates: *[]}}
     * @private
     */
    private _parseBasisLine;
    /**
     * Returns array of ALL unique elements used in basis.
     * @return {String[]}
     */
    get uniqueElements(): string[];
    /**
     * Returns array of regular bases extracted from current combinatorial basis.
     * @return {Basis[]|Object[]}
     */
    get allBasisConfigs(): Object[] | Basis[];
    /**
     * Returns array of regular bases extracted from current combinatorial basis with combinations.
     * @private
     */
    private _combination;
    /**
     * Returns array of regular bases extracted from current combinatorial basis with permutations.
     * @return {Basis[]}
     * @private
     */
    private _permutation;
    /**
     * Returns true if current combinatorial basis contains more than one regular basis.
     * @return {Boolean}
     */
    isCombinatorial(): boolean;
}
import { Basis } from "../basis/basis";
