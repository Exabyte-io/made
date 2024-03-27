"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.Basis = void 0;
// @ts-ignore
const periodic_table_js_1 = require("@exabyte-io/periodic-table.js");
const underscore_1 = __importDefault(require("underscore"));
const underscore_string_1 = __importDefault(require("underscore.string"));
const array_with_ids_1 = require("../abstract/array_with_ids");
const constants_1 = require("../constants");
const lattice_1 = require("../lattice/lattice");
const math_1 = __importDefault(require("../math"));
/**
 * A class representing a crystal basis.
 */
class Basis {
    constructor({ elements = ["Si"], coordinates = [[0, 0, 0]], units, cell = Basis.defaultCell, // by default, assume a cubic unary cell
    isEmpty = false, // whether to generate an empty Basis
    labels = [], }) {
        const _elements = isEmpty ? [] : elements;
        const _coordinates = isEmpty ? [] : coordinates;
        const _units = units || Basis.unitsOptionsDefaultValue;
        if (!units) {
            console.warn("Basis.constructor: units are not provided => set to crystal");
        }
        // assert that elements and coordinates have ids if not already passed in config + store Array helper classes
        this._elements = new array_with_ids_1.ArrayWithIds(_elements);
        this._coordinates = new array_with_ids_1.ArrayWithIds(_coordinates);
        this.units = _units;
        this.cell = cell;
        if (!underscore_1.default.isEmpty(labels)) {
            this.labels = labels;
        }
    }
    static get unitsOptionsConfig() {
        return constants_1.ATOMIC_COORD_UNITS;
    }
    static get unitsOptionsDefaultValue() {
        return constants_1.ATOMIC_COORD_UNITS.crystal;
    }
    static get defaultCell() {
        return new lattice_1.Lattice().vectorArrays;
    }
    /**
     * Serialize class instance to JSON.
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
                    6.12323399573677e-17
                ],
                [
                    1.60812264967664e-16,
                    1,
                    6.12323399573677e-17
                ],
                [
                    0,
                    0,
                    1
                ]
            ]
        }
     */
    toJSON() {
        const json = {
            elements: this.elements,
            coordinates: this.coordinates,
            units: this.units,
            cell: this.cell,
        };
        if (!underscore_1.default.isEmpty(this.labels)) {
            return JSON.parse(JSON.stringify({
                ...json,
                labels: this.labels,
            }));
        }
        return JSON.parse(JSON.stringify(json));
    }
    /**
     * Create an identical copy of the class instance.
     * @param extraContext - Extra context to be passed to the new class instance on creation.
     */
    clone(extraContext) {
        return new this.constructor({
            ...this.toJSON(),
            ...extraContext,
        });
    }
    getElementByIndex(idx) {
        return this._elements.getArrayElementByIndex(idx);
    }
    getCoordinateByIndex(idx) {
        return this._coordinates.getArrayElementByIndex(idx);
    }
    get elementsArray() {
        return this._elements.array;
    }
    get elements() {
        return this._elements.toJSON();
    }
    /**
     * Set basis elements to passed array.
     * @param elementsArray - New elements array.
     */
    set elements(elementsArray) {
        this._elements = new array_with_ids_1.ArrayWithIds(elementsArray);
    }
    getElementsAsObject() {
        return this._elements.toJSON();
    }
    get coordinates() {
        return this._coordinates.toJSON();
    }
    /**
     * Set basis elements to passed array.
     * @param {Array|ArrayWithIds} coordinatesNestedArray - New coordinates array.
     */
    set coordinates(coordinatesNestedArray) {
        this._coordinates = new array_with_ids_1.ArrayWithIds(coordinatesNestedArray);
    }
    get coordinatesAsArray() {
        return this._coordinates.array;
    }
    get isInCrystalUnits() {
        return this.units === constants_1.ATOMIC_COORD_UNITS.crystal;
    }
    get isInCartesianUnits() {
        return this.units === constants_1.ATOMIC_COORD_UNITS.cartesian;
    }
    toCartesian() {
        const unitCell = this.cell;
        if (this.units === constants_1.ATOMIC_COORD_UNITS.cartesian)
            return;
        this._coordinates.mapArrayInPlace((point) => math_1.default.multiply(point, unitCell));
        this.units = constants_1.ATOMIC_COORD_UNITS.cartesian;
    }
    toCrystal() {
        const unitCell = this.cell;
        if (this.units === constants_1.ATOMIC_COORD_UNITS.crystal)
            return;
        this._coordinates.mapArrayInPlace((point) => math_1.default.multiply(point, math_1.default.inv(unitCell)));
        this.units = constants_1.ATOMIC_COORD_UNITS.crystal;
    }
    /**
     * Asserts that all coordinates are in standardRepresentation (as explained below).
     */
    toStandardRepresentation() {
        this.toCrystal();
        this._coordinates.mapArrayInPlace((point) => point.map((x) => math_1.default.mod(x)));
    }
    /** A representation where all coordinates are within 0 and 1 in crystal units */
    get standardRepresentation() {
        const originalUnits = this.units;
        this.toStandardRepresentation();
        const result = this.toJSON();
        // preserve the original state
        if (originalUnits !== constants_1.ATOMIC_COORD_UNITS.crystal)
            this.toCartesian();
        return result;
    }
    /**
     * Add atom with a chemical element at coordinate.
     */
    addAtom({ element = "Si", coordinate = [0.5, 0.5, 0.5] }) {
        this._elements.addElement(element);
        this._coordinates.addElement(coordinate);
    }
    /**
     * Remove atom with a chemical element at coordinate either by passing the (element AND coordinate) OR id.
     */
    removeAtom({ element, coordinate, id }) {
        if (element && coordinate.length > 0) {
            this._elements.removeElement(element, id);
            this._coordinates.removeElement(coordinate, id);
        }
    }
    /**
     * Unique names (symbols) of the chemical elements basis. E.g. `['Si', 'Li']`
     */
    get uniqueElements() {
        return underscore_1.default.unique(this._elements.array);
    }
    /**
     * Returns unique chemical elements with their count sorted by electronegativity.
     * `{ "Fe": 4.0, "O": 8.0, "Li": 2.0}`.
     */
    get uniqueElementCountsSortedByElectronegativity() {
        return underscore_1.default.chain(this.elements)
            .sortBy("value")
            .sortBy((x) => (0, periodic_table_js_1.getElectronegativity)(x.value))
            .countBy((element) => element.value)
            .value();
    }
    /**
     * Returns chemical elements with their count wrt their original order in the basis.
     * Note: entries for the same element separated by another element are considered separately.
     * [{"count":1, "value":"Zr"}, {"count":23, "value":"H"}, {"count":11, "value":"Zr"}, {"count":1, "value":"H"}]
     */
    get elementCounts() {
        const elementCounts = [];
        this.getElementsAsObject().forEach((element, index) => {
            const previousElement = this.getElementsAsObject()[index - 1];
            if (previousElement && previousElement.value === element.value) {
                const previousElementCount = elementCounts[elementCounts.length - 1];
                previousElementCount.count += 1;
            }
            else {
                elementCounts.push({
                    count: 1,
                    value: element.value,
                });
            }
        });
        return elementCounts;
    }
    /**
     * Reduced formula in IUPAC format. E.g., Na2SO4
     */
    get formula() {
        const counts = this.uniqueElementCountsSortedByElectronegativity;
        const countsValues = underscore_1.default.values(counts);
        const gcd = countsValues.length > 1 ? math_1.default.gcd(...countsValues) : countsValues[0];
        return underscore_1.default.pairs(counts)
            .map((x) => x[0] + (x[1] / gcd === 1 ? "" : x[1] / gcd))
            .reduce((mem, item) => {
            return mem + item;
        }, "");
    }
    /**
     * Returns the unit cell formula as object `{ "Fe": 4.0, "O": 8.0, "Li": 2.0}`
     */
    get unitCellFormula() {
        const counts = this.uniqueElementCountsSortedByElectronegativity;
        return underscore_1.default.pairs(counts)
            .map((x) => x[0] + (x[1] === 1 ? "" : x[1]))
            .reduce((mem, item) => {
            return mem + item;
        }, "");
    }
    /**
     * Returns a nested array with elements and their corresponding coordinates
     * @example Output: [ ["Si", [0,0,0]], ["Si", [0.5,0.5,0.5]] ]
     */
    get elementsAndCoordinatesArray() {
        return this._elements.array.map((element, idx) => {
            const coordinates = this.getCoordinateByIndex(idx);
            return [element, coordinates];
        });
    }
    /**
     * Returns a nested array with elements and their corresponding coordinates with labels
     * @example Output: [ ["Si", [0,0,0], ['1']], ["Si", [0.5,0.5,0.5]] , ['2']]
     */
    get elementsAndCoordinatesAndLabelsArray() {
        return this._elements.array.map((element, idx) => {
            const coordinates = this.getCoordinateByIndex(idx);
            const atomicLabel = this.atomicLabelsArray[idx];
            return [element, coordinates, atomicLabel];
        });
    }
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
    getAsSortedString() {
        // make a copy to prevent modifying class values
        const clsInstance = new Basis(this.toJSON());
        clsInstance.toStandardRepresentation();
        const standardRep = clsInstance.elementsAndCoordinatesAndLabelsArray.map((entry) => {
            const element = entry[0];
            const coordinate = entry[1];
            const atomicLabel = entry[2];
            const toleratedCoordinates = coordinate.map((x) => math_1.default.round(x, constants_1.HASH_TOLERANCE));
            return `${element}${atomicLabel} ${toleratedCoordinates.join()}`;
        });
        return `${standardRep.sort().join(";")};`;
    }
    /**
     * Returns a string for hash calculation (in crystal units)
     */
    get hashString() {
        const originallyInCartesianUnits = this.isInCartesianUnits;
        this.toCrystal();
        const hashString = this.getAsSortedString();
        // preserve the original state
        if (originallyInCartesianUnits)
            this.toCartesian();
        return hashString;
    }
    /* Returns array of atomic labels E.g., ["1", "2", "", ""] */
    get atomicLabelsArray() {
        var _a;
        const labelsArray = Array.from({ length: this.elements.length }, (_) => "");
        (_a = this.labels) === null || _a === void 0 ? void 0 : _a.map((item) => {
            labelsArray[item.id] = item.value.toString();
        });
        return labelsArray;
    }
    /* Returns array of elements with labels E.g., ["Fe1", "Fe2", "O", "O"] */
    get elementsWithLabelsArray() {
        const elements = this.elementsArray;
        const labels = this.atomicLabelsArray;
        const elementsWithLabels = [];
        elements.forEach((symbol, idx) => elementsWithLabels.push(symbol + labels[idx]));
        return elementsWithLabels;
    }
    /**
     * Returns an array of strings with chemical elements and their atomic positions.
     * E.g., ``` ['Si 0 0 0', 'Li 0.5 0.5 0.5']```
     */
    get atomicPositions() {
        return this.elementsAndCoordinatesAndLabelsArray.map((entry, idx) => {
            const element = entry[0];
            const coordinate = entry[1];
            const atomicLabel = this.atomicLabelsArray[idx];
            return `${element}${atomicLabel} ${coordinate
                .map((x) => underscore_string_1.default.sprintf("%14.9f", x).trim())
                .join(" ")}`;
        });
    }
    /**
     * @summary Returns number of atoms in material
     */
    get nAtoms() {
        return this._elements.array.length;
    }
    // helpers
    /**
     * @summary Returns true if bases are equal, otherwise - false.
     * @param anotherBasisClsInstance {Basis} Another Basis.
     */
    isEqualTo(anotherBasisClsInstance) {
        return this.hashString === anotherBasisClsInstance.hashString;
    }
    /**
     * @summary Returns true if basis cells are equal, otherwise - false.
     * @param anotherBasisClsInstance {Basis} Another Basis.
     */
    hasEquivalentCellTo(anotherBasisClsInstance) {
        // this.cell {Array} - Cell Vectors 1, eg. [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        // prettier-ignore
        return !this.cell
            .map((vector, idx) => {
            return math_1.default.vEqualWithTolerance(vector, anotherBasisClsInstance.cell[idx]);
        })
            .some((x) => !x);
    }
    /**
     * @summary function returns the minimum basis lattice size for a structure.
     * The lattice size is based on the atomic radius of an element if the basis contains a single atom.
     * The lattice size is based on the maximum pairwise distance across a structure if the basis contains > 2 atoms.
     */
    getMinimumLatticeSize(latticeScalingFactor = lattice_1.nonPeriodicLatticeScalingFactor) {
        let latticeSizeAdditiveContribution = 0;
        if (this._elements.array.length === 1) {
            const elementSymbol = this._elements.getArrayElementByIndex(0);
            latticeSizeAdditiveContribution = (0, periodic_table_js_1.getElementAtomicRadius)(elementSymbol);
        }
        const moleculeLatticeSize = this.maxPairwiseDistance * latticeScalingFactor;
        const latticeSize = latticeSizeAdditiveContribution + moleculeLatticeSize;
        return math_1.default.precise(latticeSize, 4);
    }
    /**
     * @summary function returns an array of overlapping atoms within specified tolerance.
     */
    getOverlappingAtoms() {
        // to simplify calculations, convert to cartesian coordinates
        this.toCartesian();
        const { coordinates, elements } = this;
        const overlaps = [];
        // temporary value for overlap approximation, where atoms most certainly can't be located
        const overlapCoefficient = 0.75;
        coordinates.forEach((entry1, i) => {
            for (let j = i + 1; j < coordinates.length; j++) {
                const entry2 = coordinates[j];
                const el1 = elements[i].value;
                const el2 = elements[j].value;
                const tolerance = overlapCoefficient *
                    ((0, periodic_table_js_1.getElementAtomicRadius)(el1) + (0, periodic_table_js_1.getElementAtomicRadius)(el2)); // in angstroms
                // @ts-ignore
                const distance = math_1.default.vDist(entry1.value, entry2.value);
                if (distance < tolerance) {
                    overlaps.push({
                        id1: i,
                        id2: j,
                        element1: el1,
                        element2: el2,
                    });
                }
            }
        });
        this.toCrystal();
        return overlaps;
    }
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
    get maxPairwiseDistance() {
        const originalUnits = this.units;
        this.toCartesian();
        let maxDistance = 0;
        if (this._elements.array.length >= 2) {
            for (let i = 0; i < this._elements.array.length; i++) {
                for (let j = i + 1; j < this._elements.array.length; j++) {
                    const distance = math_1.default.vDist(this._coordinates.getArrayElementByIndex(i), this._coordinates.getArrayElementByIndex(j));
                    // @ts-ignore
                    if (distance > maxDistance) {
                        // @ts-ignore
                        maxDistance = distance;
                    }
                }
            }
        }
        if (originalUnits !== constants_1.ATOMIC_COORD_UNITS.cartesian)
            this.toCrystal();
        return math_1.default.precise(maxDistance, 4);
    }
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
    get centerOfCoordinatesPoint() {
        const transposedBasisCoordinates = math_1.default.transpose(this._coordinates.array);
        const centerOfCoordinatesVectors = [];
        for (let i = 0; i < 3; i++) {
            const center = // @ts-ignore
             transposedBasisCoordinates[i].reduce((a, b) => a + b) / this._elements.array.length;
            centerOfCoordinatesVectors.push(math_1.default.precise(center, 4));
        }
        return centerOfCoordinatesVectors;
    }
    /**
     * @summary Function translates coordinates by the vector passed as an argument.
     */
    translateByVector(translationVector) {
        // @ts-ignore
        this._coordinates.mapArrayInPlace((x) => math_1.default.add(x, translationVector));
    }
}
exports.Basis = Basis;
