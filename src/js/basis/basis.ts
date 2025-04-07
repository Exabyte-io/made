// @ts-ignore
import { getElectronegativity, getElementAtomicRadius } from "@exabyte-io/periodic-table.js";
import _ from "underscore";
import s from "underscore.string";

import { ArrayWithIds } from "../abstract/array_with_ids";
import { ObjectWithIdAndValue, ValueOrObjectArray } from "../abstract/scalar_with_id";
import { ATOMIC_COORD_UNITS, HASH_TOLERANCE } from "../constants";
import { Lattice, nonPeriodicLatticeScalingFactor } from "../lattice/lattice";
import { Vector } from "../lattice/types";
import math from "../math";
import { Coordinate } from "./types";

export interface BasisProps {
    elements: ValueOrObjectArray<string>; // chemical elements for atoms in basis.
    coordinates: ValueOrObjectArray<Coordinate>; // coordinates for the atoms in basis.
    labels?: { id: number; value: number }[];
    units: string; // units for the coordinates (eg. angstrom, crystal).
    cell: Vector[]; // crystal cell corresponding to the basis (eg. to convert to crystal coordinates).
    isEmpty?: boolean; // crystal cell corresponding to the basis (eg. to convert to crystal coordinates).
}

export interface Atom {
    id?: number; // numeric id of the element (optional).
    element: string; // Chemical element
    coordinate: Coordinate; // 3-dimensional coordinate
}

export interface ElementCount {
    count: number;
    value: string;
}

export interface BasisSchema {
    elements: ObjectWithIdAndValue<string>[];
    labels?: { id: number; value: number }[];
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
export class Basis {
    _elements: ArrayWithIds<string>;

    _coordinates: ArrayWithIds<Coordinate>;

    labels?: { id: number; value: number }[];

    units: string;

    cell: Vector[];

    constructor({
        elements = ["Si"],
        coordinates = [[0, 0, 0]],
        units,
        cell = Basis.defaultCell, // by default, assume a cubic unary cell
        isEmpty = false, // whether to generate an empty Basis
        labels = [],
    }: BasisProps) {
        const _elements = isEmpty ? [] : elements;
        const _coordinates = isEmpty ? [] : coordinates;
        const _units = units || Basis.unitsOptionsDefaultValue;

        if (!units) {
            console.warn("Basis.constructor: units are not provided => set to crystal");
        }

        // assert that elements and coordinates have ids if not already passed in config + store Array helper classes
        this._elements = new ArrayWithIds<string>(_elements);
        this._coordinates = new ArrayWithIds<Coordinate>(_coordinates);
        this.units = _units;
        this.cell = cell;

        if (!_.isEmpty(labels)) {
            this.labels = labels;
        }
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

    toJSON(skipRounding = false): BasisSchema {
        const json = {
            elements: this.elements,
            coordinates: skipRounding ? this.coordinates : this.coordinatesRounded,
            units: this.units,
        };

        if (!_.isEmpty(this.labels)) {
            return JSON.parse(
                JSON.stringify({
                    ...json,
                    labels: this.labels,
                }),
            );
        }

        return JSON.parse(JSON.stringify(json));
    }

    /** Return coordinates rounded to default precision */
    get coordinatesRounded() {
        return this.coordinates.map((coordinate) => {
            return {
                id: coordinate.id,
                value: coordinate.value.map((x) => math.precise(math.roundToZero(x))),
            };
        });
    }

    /** Return cell with vectors values rounded to default precision */
    get cellRounded() {
        return this.cell.map((vector) => vector.map((x) => math.precise(math.roundToZero(x))));
    }

    /**
     * Create an identical copy of the class instance.
     * @param extraContext - Extra context to be passed to the new class instance on creation.
     */
    clone(extraContext?: Partial<BasisProps>): Basis {
        return new (this.constructor as typeof Basis)({
            ...this.toJSON(),
            cell: this.cell,
            ...extraContext,
        });
    }

    getElementByIndex(idx: number): string {
        return this._elements.getArrayElementByIndex(idx);
    }

    getCoordinateByIndex(idx: number): Coordinate {
        return this._coordinates.getArrayElementByIndex(idx);
    }

    get elementsArray(): string[] {
        return this._elements.array;
    }

    get elements(): ObjectWithIdAndValue<string>[] {
        return this._elements.toJSON();
    }

    /**
     * Set basis elements to passed array.
     * @param elementsArray - New elements array.
     */
    set elements(elementsArray: string[] | ObjectWithIdAndValue<string>[]) {
        this._elements = new ArrayWithIds<string>(elementsArray);
    }

    getElementsAsObject(): ObjectWithIdAndValue<string>[] {
        return this._elements.toJSON();
    }

    get coordinates(): ObjectWithIdAndValue<Coordinate>[] {
        return this._coordinates.toJSON();
    }

    /**
     * Set basis elements to passed array.
     * @param {Array|ArrayWithIds} coordinatesNestedArray - New coordinates array.
     */
    set coordinates(coordinatesNestedArray: Coordinate[] | ObjectWithIdAndValue<Coordinate>[]) {
        this._coordinates = new ArrayWithIds<Coordinate>(coordinatesNestedArray);
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
        this._coordinates.mapArrayInPlace(
            (point) => math.multiply(point, unitCell) as unknown as Coordinate,
        );
        this.units = ATOMIC_COORD_UNITS.cartesian;
    }

    toCrystal() {
        const unitCell = this.cell;
        if (this.units === ATOMIC_COORD_UNITS.crystal) return;
        this._coordinates.mapArrayInPlace(
            (point) => math.multiply(point, math.inv(unitCell)) as unknown as Coordinate,
        );
        this.units = ATOMIC_COORD_UNITS.crystal;
    }

    /**
     * Asserts that all coordinates are in standardRepresentation (as explained below).
     */
    toStandardRepresentation() {
        this.toCrystal();
        this._coordinates.mapArrayInPlace((point) => point.map((x) => math.mod(x)) as Coordinate);
    }

    /** A representation where all coordinates are within 0 and 1 in crystal units */
    get standardRepresentation(): BasisSchema {
        const originalUnits = this.units;

        this.toStandardRepresentation();
        const result = this.toJSON();

        // preserve the original state
        if (originalUnits !== ATOMIC_COORD_UNITS.crystal) this.toCartesian();

        return result;
    }

    /**
     * Add atom with a chemical element at coordinate.
     */
    addAtom({ element = "Si", coordinate = [0.5, 0.5, 0.5] }: Atom) {
        this._elements.addElement(element);
        this._coordinates.addElement(coordinate);
    }

    /**
     * Remove atom with a chemical element at coordinate either by passing the (element AND coordinate) OR id.
     */
    removeAtom({ element, coordinate, id }: Atom) {
        if (element && coordinate.length > 0) {
            this._elements.removeElement(element, id);
            this._coordinates.removeElement(coordinate, id);
        }
    }

    /**
     * Unique names (symbols) of the chemical elements basis. E.g. `['Si', 'Li']`
     */
    get uniqueElements(): string[] {
        return _.unique(this._elements.array);
    }

    /**
     * Returns unique chemical elements with their count sorted by electronegativity.
     * `{ "Fe": 4.0, "O": 8.0, "Li": 2.0}`.
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
     */
    get elementCounts(): ElementCount[] {
        const elementCounts: ElementCount[] = [];
        this.getElementsAsObject().forEach((element, index) => {
            const previousElement = this.getElementsAsObject()[index - 1];
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
     */
    get formula(): string {
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
     */
    get unitCellFormula(): string {
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
     */
    get elementsAndCoordinatesArray(): [string, Coordinate][] {
        return this._elements.array.map((element, idx) => {
            const coordinates = this.getCoordinateByIndex(idx);
            return [element, coordinates];
        });
    }

    /**
     * Returns a nested array with elements and their corresponding coordinates with labels
     * @example Output: [ ["Si", [0,0,0], ['1']], ["Si", [0.5,0.5,0.5]] , ['2']]
     */
    get elementsAndCoordinatesAndLabelsArray(): [string, Coordinate, string][] {
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
    getAsSortedString(): string {
        // make a copy to prevent modifying class values
        const clsInstance = new Basis(this.toJSON());
        clsInstance.toStandardRepresentation();
        const standardRep = clsInstance.elementsAndCoordinatesAndLabelsArray.map((entry) => {
            const element = entry[0];
            const coordinate = entry[1];
            const atomicLabel = entry[2];
            const toleratedCoordinates = coordinate.map((x) => math.round(x, HASH_TOLERANCE));
            return `${element}${atomicLabel} ${toleratedCoordinates.join()}`;
        });
        return `${standardRep.sort().join(";")};`;
    }

    /**
     * Returns a string for hash calculation (in crystal units)
     */
    get hashString(): string {
        const originallyInCartesianUnits = this.isInCartesianUnits;

        this.toCrystal();
        const hashString = this.getAsSortedString();

        // preserve the original state
        if (originallyInCartesianUnits) this.toCartesian();

        return hashString;
    }

    /* Returns array of atomic labels E.g., ["1", "2", "", ""] */
    get atomicLabelsArray(): string[] {
        const labelsArray = Array.from({ length: this.elements.length }, (_) => "");
        // https://dev.to/maafaishal/benchmarking-for-while-forof-and-arrayforeach-using-performancenow-1jjg
        if (this.labels?.length) {
            for (let i = 0; i < this.labels.length; i++) {
                labelsArray[this.labels[i].id] = this.labels[i].value.toString();
            }
        }
        return labelsArray;
    }

    /* Returns array of elements with labels E.g., ["Fe1", "Fe2", "O", "O"] */
    get elementsWithLabelsArray(): string[] {
        const elements = this.elementsArray;
        const labels = this.atomicLabelsArray;
        const elementsWithLabels: string[] = [];
        elements.forEach((symbol, idx) => elementsWithLabels.push(symbol + labels[idx]));
        return elementsWithLabels;
    }

    /**
     * Returns an array of strings with chemical elements and their atomic positions.
     * E.g., ``` ['Si 0 0 0', 'Li 0.5 0.5 0.5']```
     */
    get atomicPositions(): string[] {
        return this.elementsAndCoordinatesAndLabelsArray.map((entry, idx) => {
            const element = entry[0];
            const coordinate = entry[1];
            const atomicLabel = this.atomicLabelsArray[idx];
            return `${element}${atomicLabel} ${coordinate
                .map((x) => s.sprintf("%14.9f", x).trim())
                .join(" ")}`;
        });
    }

    /**
     * @summary Returns number of atoms in material
     */
    get nAtoms(): number {
        return this._elements.array.length;
    }

    // helpers

    /**
     * @summary Returns true if bases are equal, otherwise - false.
     * @param anotherBasisClsInstance {Basis} Another Basis.
     */
    isEqualTo(anotherBasisClsInstance: Basis): boolean {
        return this.hashString === anotherBasisClsInstance.hashString;
    }

    /**
     * @summary Returns true if basis cells are equal, otherwise - false.
     * @param anotherBasisClsInstance {Basis} Another Basis.
     */
    hasEquivalentCellTo(anotherBasisClsInstance: Basis): boolean {
        // this.cell {Array} - Cell Vectors 1, eg. [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        // prettier-ignore
        return !this.cell
            .map((vector, idx) => {
                return math.vEqualWithTolerance(vector, anotherBasisClsInstance.cell[idx]);
            })
            .some((x) => !x);
    }

    /**
     * @summary function returns the minimum basis lattice size for a structure.
     * The lattice size is based on the atomic radius of an element if the basis contains a single atom.
     * The lattice size is based on the maximum pairwise distance across a structure if the basis contains > 2 atoms.
     */
    getMinimumLatticeSize(latticeScalingFactor = nonPeriodicLatticeScalingFactor): number {
        let latticeSizeAdditiveContribution = 0;
        if (this._elements.array.length === 1) {
            const elementSymbol = this._elements.getArrayElementByIndex(0);
            latticeSizeAdditiveContribution = getElementAtomicRadius(elementSymbol);
        }
        const moleculeLatticeSize = this.maxPairwiseDistance * latticeScalingFactor;
        const latticeSize = latticeSizeAdditiveContribution + moleculeLatticeSize;
        return math.precise(latticeSize, 4);
    }

    /**
     * @summary function returns an array of overlapping atoms within specified tolerance.
     */
    getOverlappingAtoms(): Overlap[] {
        // to simplify calculations, convert to cartesian coordinates
        this.toCartesian();
        const { coordinates, elements } = this;
        const overlaps: Overlap[] = [];
        // temporary value for overlap approximation, where atoms most certainly can't be located
        const overlapCoefficient = 0.75;

        coordinates.forEach((entry1, i) => {
            for (let j = i + 1; j < coordinates.length; j++) {
                const entry2 = coordinates[j];
                const el1 = elements[i].value;
                const el2 = elements[j].value;

                const tolerance =
                    overlapCoefficient *
                    (getElementAtomicRadius(el1) + getElementAtomicRadius(el2)); // in angstroms

                // @ts-ignore
                const distance = math.vDist(entry1.value, entry2.value) as number;
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
    get maxPairwiseDistance(): number {
        const originalUnits = this.units;
        this.toCartesian();
        let maxDistance = 0;
        if (this._elements.array.length >= 2) {
            for (let i = 0; i < this._elements.array.length; i++) {
                for (let j = i + 1; j < this._elements.array.length; j++) {
                    const distance = math.vDist(
                        this._coordinates.getArrayElementByIndex(i),
                        this._coordinates.getArrayElementByIndex(j),
                    );
                    // @ts-ignore
                    if (distance > maxDistance) {
                        // @ts-ignore
                        maxDistance = distance;
                    }
                }
            }
        }
        if (originalUnits !== ATOMIC_COORD_UNITS.cartesian) this.toCrystal();
        return math.precise(maxDistance, 4);
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
    get centerOfCoordinatesPoint(): number[] {
        const transposedBasisCoordinates = math.transpose(this._coordinates.array);
        const centerOfCoordinatesVectors = [];
        for (let i = 0; i < 3; i++) {
            const center = // @ts-ignore
                transposedBasisCoordinates[i].reduce((a, b) => a + b) / this._elements.array.length;
            centerOfCoordinatesVectors.push(math.precise(center, 4));
        }
        return centerOfCoordinatesVectors;
    }

    /**
     * @summary Function translates coordinates by the vector passed as an argument.
     */
    translateByVector(translationVector: number[]) {
        // @ts-ignore
        this._coordinates.mapArrayInPlace((x) => math.add(x, translationVector));
    }
}
