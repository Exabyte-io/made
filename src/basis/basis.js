import _ from "underscore";
import s from "underscore.string";
import {getElectronegativity} from "periodic-table";

import math from "../math";
import {ArrayWithIds} from "../primitive";
import {Lattice} from "../lattice/lattice";
import {ATOMIC_COORD_UNITS, HASH_TOLERANCE} from "../constants";

export class Basis {
    constructor({
                    elements = ["Si"],
                    coordinates = [[0, 0, 0]],
                    units,
                    cell = Basis.defaultCell,   // by default assume a cubic unary cell
                    isEmpty = false,            // whether to generate an empty Basis
                }) {
        if (!units) {
            units = Basis.unitsOptionsDefaultValue;
            console.warn("Basis.constructor: units are not provided => set to crystal");
        }
        if (isEmpty) {
            elements = coordinates = [];
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

    set elements(elementsArray) {
        this._elements = new ArrayWithIds(elementsArray);
    }

    get coordinates() {
        return this._coordinates.toJSON();
    }

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
        this._coordinates.mapArrayInPlace(point => math.multiply(point, unitCell));
        this.units = ATOMIC_COORD_UNITS.cartesian;
    }

    toCrystal() {
        const unitCell = this.cell;
        if (this.units === ATOMIC_COORD_UNITS.crystal) return;
        this._coordinates.mapArrayInPlace(point => math.multiply(point, math.inv(unitCell)));
        this.units = ATOMIC_COORD_UNITS.crystal;
    }

    toJSON() {
        return JSON.parse(JSON.stringify(_.pick(this, ["elements", "coordinates", "units", "cell"])));
    }

    clone(extraContext) {
        return new this.constructor(Object.assign({}, this.toJSON(), extraContext));
    }

    // assert that all coordinates are within 0 and 1 in crystal units
    toStandardRepresentation() {
        this.toCrystal();
        this._coordinates.mapArrayInPlace(point => point.map(x => math.mod(x)));
    }

    get standardRepresentation() {
        const originalUnits = this.units;

        this.toStandardRepresentation();
        const result = this.toJSON();

        // preserve the original state
        if (originalUnits !== ATOMIC_COORD_UNITS.crystal) this.toCartesian();

        return result;
    }

    addAtom({element = "Si", coordinate = [0.5, 0.5, 0.5]}) {
        this._elements.addElement(element);
        this._coordinates.addElement(coordinate);
    }

    removeAtom({element, coordinate, id}) {
        if (element && coordinate.length > 0) {
            this._elements.removeElement(element, id);
            this._coordinates.removeElement(coordinate, id);
        }
    }

    removeAtomAtCoordinate({coordinate, id}) {
        if (coordinate.length > 0) {
            const index = this._coordinates.getArrayIndexByPredicate(
                arrayCoordinate => (math.roundToZero(math.vDist(arrayCoordinate, coordinate)) === 0)
            );
            (index > -1) && this._elements.removeElement(null, index);
            (index > -1) && this._coordinates.removeElement(null, index);
        }
    }

    /**
     * @summary Unique element names used in the material. E.g. `['Si', 'Li']`
     */
    get uniqueElements() {
        return _.unique(this._elements.array);
    }

    /**
     * @summary Returns unit cell formula as object `{ "Fe": 4.0, "O": 8.0, "Li": 2.0}`
     */
    get elementCounts() {
        return _.chain(this.elements)
            .sortBy('value').sortBy(x => getElectronegativity(x.value))
            .countBy(element => element.value).value();
    }

    /**
     * @summary Returns a reduced formula in IUPAC format. E.g., Na2SO4
     * TODO: sort elements by electronegativity
     */
    get formula() {
        const counts = this.elementCounts;
        const countsValues = _.values(counts);
        const gcd = countsValues.length > 1 ? math.gcd.apply(math, countsValues) : countsValues[0];
        return _.pairs(counts).map(x => x[0] + ((x[1] / gcd) === 1 ? '' : (x[1] / gcd))).reduce((mem, item) => {
            return mem + item;
        }, '');
    }

    /**
     * @summary Returns unit cell formula as object `{ "Fe": 4.0, "O": 8.0, "Li": 2.0}`
     */
    get unitCellFormula() {
        const counts = this.elementCounts;
        return _.pairs(counts).map(x => x[0] + (x[1] === 1 ? '' : x[1])).reduce((mem, item) => {
            return mem + item;
        }, '');
    }

    /**
     * @summary Returns a nested array with elements and their corresponding coordinates
     * @example Output: [ ["Si", [0,0,0]], ["Si", [0.5,0.5,0.5]] ]
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
     *     result: "Si 0,1,0;Si 1,0,0"
     * ```
     * This guarantees independence of hash from elements array ordering.
     * @return {String}
     */
    getAsSortedString(inStandardRepresentation = true) {
        // make a copy to prevent modifying class values
        const clsInstance = new Basis(this.toJSON());
        clsInstance.toStandardRepresentation();
        return clsInstance.elementsAndCoordinatesArray.map((entry) => {
            const element = entry[0];
            const coordinate = entry[1];
            const toleratedCoordinates = coordinate.map(x => math.round(x, HASH_TOLERANCE));
            return `${element} ${toleratedCoordinates.join()}`;
        }).sort().join(';') + ";";
    }

    // returns a string for hash calculation (in crystal units)
    get hashString() {
        const originallyInCartesianUnits = this.isInCartesianUnits;

        this.toCrystal();
        const hashString = this.getAsSortedString();

        // preserve the original state
        if (originallyInCartesianUnits) this.toCartesian();

        return hashString;
    }

    // helper to compare basis instances
    isEqualTo(anotherBasisCls) {
        return this.hashString === anotherBasisCls.hashString;
    }

    /**
     * @summary Returns atomic positions array.
     * E.g., ``` ['Si 0 0 0', 'Li 0.5 0.5 0.5']```
     * @return {String[]}
     */
    get atomicPositions() {
        const clsInstance = this;
        return this.elementsAndCoordinatesArray.map((entry) => {
            const element = entry[0];
            const coordinate = entry[1];
            return element + ' ' + coordinate.map(x => s.sprintf('%14.9f', x).trim()).join(' ');
        });
    }

    /**
     * @summary Returns number of atoms in material
     */
    get nAtoms() {
        return this._elements.array.length;
    }

}
