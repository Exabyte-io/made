import _ from "underscore";
import { ArrayWithIds } from "../abstract/array_with_ids";
import { ObjectWithIdAndValue, ValueOrObjectArray } from "../abstract/scalar_with_id";
import { ATOMIC_COORD_UNITS } from "../constants";
import { Vector } from "../lattice/types";
import { Coordinate } from "./types";
export interface BasisProps {
    elements: ValueOrObjectArray<string>;
    coordinates: ValueOrObjectArray<Coordinate>;
    units: string;
    cell: Vector[];
    isEmpty?: boolean;
}
export interface Atom {
    id?: number;
    element: string;
    coordinate: Coordinate;
}
export interface ElementCount {
    count: number;
    value: string;
}
export interface BasisSchema {
    elements: ObjectWithIdAndValue<string>[];
    coordinates: ObjectWithIdAndValue<Coordinate>[];
    units: string;
    cell: Vector[];
}
interface Overlap {
    id1: number;
    id2: number;
    element1: string;
    element2: string;
}
/**
 * A class representing a crystal basis.
 */
export declare class Basis {
    _elements: ArrayWithIds<string>;
    _coordinates: ArrayWithIds<Coordinate>;
    units: string;
    cell: Vector[];
    precision: number;
    constructor({ elements, coordinates, units, cell, // by default, assume a cubic unary cell
    isEmpty, }: BasisProps);
    static get unitsOptionsConfig(): typeof ATOMIC_COORD_UNITS;
    static get unitsOptionsDefaultValue(): string;
    static get defaultCell(): [import("@mat3ra/esse/lib/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/lib/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/lib/js/types").ArrayOf3NumberElementsSchema];
    /**
     * Serialize class instance to JSON.
     * @param skipRounding - Whether to skip rounding the resulting lattice values, defaults to `false`.
     * @example As below:
     {
            "elements" : [
                {
                    "id" : 0,
                    "value" : "Si"
                },
                {
                    "id" : 1,
                    "value" : "Si"
                }
            ],
            "coordinates" : [
                {
                    "id" : 0,
                    "value" : [
                        0,
                        0,
                        0
                    ]
                },
                {
                    "id" : 1,
                    "value" : [
                        0.25,
                        0.25,
                        0.25
                    ]
                }
            ],
            "units" : "crystal",
            "cell" : [
                [
                    1,
                    0,
                    0
                ],
                [
                    0,
                    1,
                    0
                ],
                [
                    0,
                    0,
                    1
                ]
            ]
        }
     */
    toJSON(skipRounding?: boolean): BasisSchema;
    /** Round coordinates to the specified precision */
    roundCoordinates(): void;
    /** Round cell values to the specified precision */
    roundCell(): void;
    /**
     * Create an identical copy of the class instance.
     * @param extraContext - Extra context to be passed to the new class instance on creation.
     */
    clone(extraContext?: Partial<BasisProps>): Basis;
    getElementByIndex(idx: number): string;
    getCoordinateByIndex(idx: number): Coordinate;
    get elementsArray(): string[];
    get elements(): ObjectWithIdAndValue<string>[];
    /**
     * Set basis elements to passed array.
     * @param elementsArray - New elements array.
     */
    set elements(elementsArray: string[] | ObjectWithIdAndValue<string>[]);
    getElementsAsObject(): ObjectWithIdAndValue<string>[];
    get coordinates(): ObjectWithIdAndValue<Coordinate>[];
    /**
     * Set basis elements to passed array.
     * @param {Array|ArrayWithIds} coordinatesNestedArray - New coordinates array.
     */
    set coordinates(coordinatesNestedArray: Coordinate[] | ObjectWithIdAndValue<Coordinate>[]);
    get coordinatesAsArray(): Coordinate[];
    get isInCrystalUnits(): boolean;
    get isInCartesianUnits(): boolean;
    toCartesian(): void;
    toCrystal(): void;
    /**
     * Asserts that all coordinates are in standardRepresentation (as explained below).
     */
    toStandardRepresentation(): void;
    /** A representation where all coordinates are within 0 and 1 in crystal units */
    get standardRepresentation(): BasisSchema;
    /**
     * Add atom with a chemical element at coordinate.
     */
    addAtom({ element, coordinate }: Atom): void;
    /**
     * Remove atom with a chemical element at coordinate either by passing the (element AND coordinate) OR id.
     */
    removeAtom({ element, coordinate, id }: Atom): void;
    /**
     * Unique names (symbols) of the chemical elements basis. E.g. `['Si', 'Li']`
     */
    get uniqueElements(): string[];
    /**
     * Returns unique chemical elements with their count sorted by electronegativity.
     * `{ "Fe": 4.0, "O": 8.0, "Li": 2.0}`.
     */
    get uniqueElementCountsSortedByElectronegativity(): _.Dictionary<number>;
    /**
     * Returns chemical elements with their count wrt their original order in the basis.
     * Note: entries for the same element separated by another element are considered separately.
     * [{"count":1, "value":"Zr"}, {"count":23, "value":"H"}, {"count":11, "value":"Zr"}, {"count":1, "value":"H"}]
     */
    get elementCounts(): ElementCount[];
    /**
     * Reduced formula in IUPAC format. E.g., Na2SO4
     */
    get formula(): string;
    /**
     * Returns the unit cell formula as object `{ "Fe": 4.0, "O": 8.0, "Li": 2.0}`
     */
    get unitCellFormula(): string;
    /**
     * Returns a nested array with elements and their corresponding coordinates
     * @example Output: [ ["Si", [0,0,0]], ["Si", [0.5,0.5,0.5]] ]
     */
    get elementsAndCoordinatesArray(): [string, Coordinate][];
    /**
     * @summary Concatenates elements and sorts them in alphanumeric order.
     * E.g.,
     * ```
     *     elements: [{id: 0, value: 'Si'}, {id: 1, value: 'Si'}]
     *     coordinates: [{id: 0, value: [1,0,0]}, {id: 1, value: [0, 1, 0]}]
     *
     *     result: "Si 0,1,0;Si 1,0,0"
     * ```
     * This guarantees the independence from the order in the elements array.
     */
    getAsSortedString(): string;
    /**
     * Returns a string for hash calculation (in crystal units)
     */
    get hashString(): string;
    /**
     * Returns an array of strings with chemical elements and their atomic positions.
     * E.g., ``` ['Si 0 0 0', 'Li 0.5 0.5 0.5']```
     */
    get atomicPositions(): string[];
    /**
     * @summary Returns number of atoms in material
     */
    get nAtoms(): number;
    /**
     * @summary Returns true if bases are equal, otherwise - false.
     * @param anotherBasisClsInstance {Basis} Another Basis.
     */
    isEqualTo(anotherBasisClsInstance: Basis): boolean;
    /**
     * @summary Returns true if basis cells are equal, otherwise - false.
     * @param anotherBasisClsInstance {Basis} Another Basis.
     */
    hasEquivalentCellTo(anotherBasisClsInstance: Basis): boolean;
    /**
     * @summary function returns the minimum basis lattice size for a structure.
     * The lattice size is based on the atomic radius of an element if the basis contains a single atom.
     * The lattice size is based on the maximum pairwise distance across a structure if the basis contains > 2 atoms.
     */
    getMinimumLatticeSize(latticeScalingFactor?: number): number;
    /**
     * @summary function returns an array of overlapping atoms within specified tolerance.
     */
    getOverlappingAtoms(): Overlap[];
    /**
     * @summary function returns the max distance between pairs of basis coordinates by
     * calculating the distance between pairs of basis coordinates.
     * basis coordinates = [[x1, y1, z1], [x2, y2, z2], ... [xn, yn, zn]]
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
     */
    get maxPairwiseDistance(): number;
    /**
     * @summary Function takes basis coordinates and transposes them so that the values for each dimension of the
     *  the basis are in their own nested array.
     *  Then the center point for each dimension of the coordinates is calculated.
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
     */
    get centerOfCoordinatesPoint(): number[];
    /**
     * @summary Function translates coordinates by the vector passed as an argument.
     */
    translateByVector(translationVector: number[]): void;
}
export {};
