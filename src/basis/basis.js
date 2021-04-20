import { getElectronegativity } from "@exabyte-io/periodic-table.js";
import _ from "underscore";
import s from "underscore.string";

import { ArrayWithIds } from "../abstract/array_with_ids";
import { ATOMIC_COORD_UNITS, HASH_TOLERANCE } from "../constants";
import { Lattice } from "../lattice/lattice";
import math from "../math";

/**
 * A class representing a crystal basis.
 */
export class Basis {
    /**
     * Create a Basis class.
     * @param {Object} - Config object.
     * @param {Array|ArrayWithIds} elements - chemical elements for atoms in basis.
     * @param {Array|ArrayWithIds} coordinates - coordinates for the atoms in basis.
     * @param {String} units - units for the coordinates (eg. angstrom, crystal).
     * @param {Cell} cell - crystal cell corresponding to the basis (eg. to convert to crystal coordinates).
     * @param {Boolean} isEmpty - crystal cell corresponding to the basis (eg. to convert to crystal coordinates).
     */
    constructor({
        elements = ["Si"],
        coordinates = [[0, 0, 0]],
        units,
        cell = Basis.defaultCell, // by default, assume a cubic unary cell
        isEmpty = false, // whether to generate an empty Basis
    }) {
        if (!units) {
            units = Basis.unitsOptionsDefaultValue;
            console.warn("Basis.constructor: units are not provided => set to crystal");
        }
        if (isEmpty) {
            elements = [];
            coordinates = [];
        }
        // assert that elements and coordinates have ids if not already passed in config + store Array helper classes
        this._elements = new ArrayWithIds(elements);
        this._coordinates = new ArrayWithIds(coordinates);
        this.units = units;
        this.cell = cell;
    }

    static get unitsOptionsConfig() {
        return ATOMIC_COORD_UNITS;
    }

    static get unitsOptionsDefaultValue() {
        return ATOMIC_COORD_UNITS.crystal;
    }

    static get defaultCell() {
        return new Lattice().vectorArrays;
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
        return JSON.parse(
            JSON.stringify(_.pick(this, ["elements", "coordinates", "units", "cell"])),
        );
    }

    /**
     * Create an identical copy of the class instance.
     * @param {Object} extraContext - Extra context to be passed to the new class instance on creation.
     */
    clone(extraContext) {
        return new this.constructor({ ...this.toJSON(), ...extraContext });
    }

    getElementByIndex(idx) {
        return this._elements.getArrayElementByIndex(idx);
    }

    getCoordinateByIndex(idx) {
        return this._coordinates.getArrayElementByIndex(idx);
    }

    get elements() {
        return this._elements.toJSON();
    }

    get elementsArray() {
        return this._elements.array;
    }

    /**
     * Set basis elements to passed array.
     * @param {Array|ArrayWithIds} elementsArray - New elements array.
     */
    set elements(elementsArray) {
        this._elements = new ArrayWithIds(elementsArray);
    }

    get coordinates() {
        return this._coordinates.toJSON();
    }

    /**
     * Set basis elements to passed array.
     * @param {Array|ArrayWithIds} coordinatesNestedArray - New coordinates array.
     */
    set coordinates(coordinatesNestedArray) {
        this._coordinates = new ArrayWithIds(coordinatesNestedArray);
    }

    get coordinatesAsArray() {
        return this._coordinates.array;
    }

    get isInCrystalUnits() {
        return this.units === ATOMIC_COORD_UNITS.crystal;
    }

    get isInCartesianUnits() {
        return this.units === ATOMIC_COORD_UNITS.cartesian;
    }

    toCartesian() {
        const unitCell = this.cell;
        if (this.units === ATOMIC_COORD_UNITS.cartesian) return;
        this._coordinates.mapArrayInPlace((point) => math.multiply(point, unitCell));
        this.units = ATOMIC_COORD_UNITS.cartesian;
    }

    toCrystal() {
        const unitCell = this.cell;
        if (this.units === ATOMIC_COORD_UNITS.crystal) return;
        this._coordinates.mapArrayInPlace((point) => math.multiply(point, math.inv(unitCell)));
        this.units = ATOMIC_COORD_UNITS.crystal;
    }

    /**
     * Asserts that all coordinates are in standardRepresentation (as explained below).
     */
    toStandardRepresentation() {
        this.toCrystal();
        this._coordinates.mapArrayInPlace((point) => point.map((x) => math.mod(x)));
    }

    /** A representation where all coordinates are within 0 and 1 in crystal units */
    get standardRepresentation() {
        const originalUnits = this.units;

        this.toStandardRepresentation();
        const result = this.toJSON();

        // preserve the original state
        if (originalUnits !== ATOMIC_COORD_UNITS.crystal) this.toCartesian();

        return result;
    }

    /**
     * Add atom with a chemical element at coordinate.
     * @param {Object} config
     * @param {String} element - Chemical element.
     * @param {Array} coordinate - 3-dimensional coordinate.
     */
    addAtom({ element = "Si", coordinate = [0.5, 0.5, 0.5] }) {
        this._elements.addElement(element);
        this._coordinates.addElement(coordinate);
    }

    /**
     * Remove atom with a chemical element at coordinate either by passing the (element AND coordinate) OR id.
     * @param {Object} config
     * @param {String} element - Chemical element.
     * @param {Array} coordinate - 3-dimensional coordinate.
     * @param {Number} id - numeric id of the element (optional).
     */
    removeAtom({ element, coordinate, id }) {
        if (element && coordinate.length > 0) {
            this._elements.removeElement(element, id);
            this._coordinates.removeElement(coordinate, id);
        }
    }

    /**
     * Remove atom at a particular coordinate.
     * @param {Object} config
     * @param {Array} coordinate - 3-dimensional coordinate.
     */
    removeAtomAtCoordinate({ coordinate }) {
        if (coordinate.length > 0) {
            const index = this._coordinates.getArrayIndexByPredicate(
                (arrayCoordinate) =>
                    math.roundToZero(math.vDist(arrayCoordinate, coordinate)) === 0,
            );
            index > -1 && this._elements.removeElement(null, index);
            index > -1 && this._coordinates.removeElement(null, index);
        }
    }

    /**
     * Unique names (symbols) of the chemical elements basis. E.g. `['Si', 'Li']`
     * @return {Array}
     */
    get uniqueElements() {
        return _.unique(this._elements.array);
    }

    /**
     * Returns unique chemical elements with their count sorted by electronegativity.
     * `{ "Fe": 4.0, "O": 8.0, "Li": 2.0}`.
     * @return {Object}
     */
    get uniqueElementCountsSortedByElectronegativity() {
        return _.chain(this.elements)
            .sortBy("value")
            .sortBy((x) => getElectronegativity(x.value))
            .countBy((element) => element.value)
            .value();
    }

    /**
     * Returns chemical elements with their count wrt their original order in the basis.
     * Note: entries for the same element separated by another element are considered separately.
     * [{"count":1, "value":"Zr"}, {"count":23, "value":"H"}, {"count":11, "value":"Zr"}, {"count":1, "value":"H"}]
     * @return {Object}
     */
    get elementCounts() {
        const elementCounts = [];
        this.elements.forEach((element, index) => {
            const previousElement = this.elements[index - 1];
            if (previousElement && previousElement.value === element.value) {
                const previousElementCount = elementCounts[elementCounts.length - 1];
                previousElementCount.count += 1;
            } else {
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
     * @return {String}
     */
    get formula() {
        const counts = this.uniqueElementCountsSortedByElectronegativity;
        const countsValues = _.values(counts);
        const gcd = countsValues.length > 1 ? math.gcd(...countsValues) : countsValues[0];
        return _.pairs(counts)
            .map((x) => x[0] + (x[1] / gcd === 1 ? "" : x[1] / gcd))
            .reduce((mem, item) => {
                return mem + item;
            }, "");
    }

    /**
     * Returns the unit cell formula as object `{ "Fe": 4.0, "O": 8.0, "Li": 2.0}`
     * @return {Array}
     */
    get unitCellFormula() {
        const counts = this.uniqueElementCountsSortedByElectronegativity;
        return _.pairs(counts)
            .map((x) => x[0] + (x[1] === 1 ? "" : x[1]))
            .reduce((mem, item) => {
                return mem + item;
            }, "");
    }

    /**
     * Returns a nested array with elements and their corresponding coordinates
     * @example Output: [ ["Si", [0,0,0]], ["Si", [0.5,0.5,0.5]] ]
     * @return {Array}
     */
    get elementsAndCoordinatesArray() {
        const clsInstance = this;
        return this._elements.array.map((element, idx) => {
            const coordinates = clsInstance.getCoordinateByIndex(idx);
            return [element, coordinates];
        });
    }

    /**
     * @summary Concatenates elements and sorts them in alphanumeric order.
     * E.g.,
     * ```
     *     elements: [{id: 1, value: 'Si'}, {id: 2, value: 'Si'}]
     *     coordinates: [{id: 1, value: [1,0,0]}, {id: 2, value: [0, 1, 0]}]
     *
     *     result: "Si 0,1,0;Si 1,0,0"
     * ```
     * This guarantees the independence from the order in the elements array.
     * @return {String}
     */
    getAsSortedString() {
        // make a copy to prevent modifying class values
        const clsInstance = new Basis(this.toJSON());
        clsInstance.toStandardRepresentation();
        return `${clsInstance.elementsAndCoordinatesArray
            .map((entry) => {
                const element = entry[0];
                const coordinate = entry[1];
                const toleratedCoordinates = coordinate.map((x) => math.round(x, HASH_TOLERANCE));
                return `${element} ${toleratedCoordinates.join()}`;
            })
            .sort()
            .join(";")};`;
    }

    /**
     * Returns a string for hash calculation (in crystal units)
     * @return {String}
     */
    get hashString() {
        const originallyInCartesianUnits = this.isInCartesianUnits;

        this.toCrystal();
        const hashString = this.getAsSortedString();

        // preserve the original state
        if (originallyInCartesianUnits) this.toCartesian();

        return hashString;
    }

    /**
     * Returns an array of strings with chemical elements and their atomic positions.
     * E.g., ``` ['Si 0 0 0', 'Li 0.5 0.5 0.5']```
     * @return {String[]}
     */
    get atomicPositions() {
        return this.elementsAndCoordinatesArray.map((entry) => {
            const element = entry[0];
            const coordinate = entry[1];
            return `${element} ${coordinate.map((x) => s.sprintf("%14.9f", x).trim()).join(" ")}`;
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
     * @return {Boolean}
     */
    isEqualTo(anotherBasisClsInstance) {
        return this.hashString === anotherBasisClsInstance.hashString;
    }

    /**
     * @summary Returns true if basis cells are equal, otherwise - false.
     * @param anotherBasisClsInstance {Basis} Another Basis.
     * @return {Boolean}
     */
    hasEquivalentCellTo(anotherBasisClsInstance) {
        // this.cell {Array} - Cell Vectors 1, eg. [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        return !this.cell
            .map((vector, idx) =>
                math.vEqualWithTolerance(vector, anotherBasisClsInstance.cell[idx]),
            )
            .some((x) => !x);
    }
}
