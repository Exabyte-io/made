// @ts-ignore
import { getElectronegativity, getElementAtomicRadius } from "@exabyte-io/periodic-table.js";
import { ArrayWithIds } from "@mat3ra/code";
import { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import { AnyObject } from "@mat3ra/esse/dist/js/esse/types";
import { BasisSchema, PointSchema } from "@mat3ra/esse/dist/js/types";
import * as _ from "underscore";
import * as s from "underscore.string";

import { Cell } from "../cell/cell";
import { ATOMIC_COORD_UNITS, HASH_TOLERANCE } from "../constants";
import { nonPeriodicLatticeScalingFactor } from "../lattice/lattice";
import math from "../math";
import { Coordinate, Coordinates } from "./coordinates";
import { Elements } from "./elements";

export interface Atom {
    id?: number; // numeric id of the element (optional).
    element: string; // Chemical element
    coordinate: Coordinate; // Coordinates of the element
}

export interface ElementCount {
    count: number;
    value: string;
}

interface Overlap {
    id1: number;
    id2: number;
    element1: string;
    element2: string;
}

export interface BasisConfig extends BasisSchema {
    cell?: Cell;
    isEmpty?: boolean;
}

const DEFAULT_BASIS_CONFIG = {
    elements: [
        {
            id: 0,
            value: "Si",
        },
        {
            id: 1,
            value: "Si",
        },
    ],
    coordinates: [
        {
            id: 0,
            value: [0, 0, 0],
        },
        {
            id: 1,
            value: [0.25, 0.25, 0.25],
        },
    ],
    units: "crystal",
};

export class Basis extends InMemoryEntity implements BasisSchema {
    static defaultConfig: BasisSchema = DEFAULT_BASIS_CONFIG as BasisSchema;

    units: BasisSchema["units"];

    cell: Cell;

    elements: BasisSchema["elements"];

    coordinates: BasisSchema["coordinates"];

    labels: BasisSchema["labels"];

    _elements: Elements;

    _coordinates: Coordinates;

    _labels: ArrayWithIds<string>;

    static fromElementsAndCoordinates({
        elements = [],
        coordinates = [],
        units = ATOMIC_COORD_UNITS.crystal as BasisSchema["units"],
        cell = new Cell(),
        labels = [],
    }): Basis {
        const elementsArrayWithIdsJSON = Elements.fromValues(elements).toJSON();
        const coordinatesArrayWithIdsJSON = Coordinates.fromValues(coordinates).toJSON();
        const labelsArrayWithIdsJSON = ArrayWithIds.fromValues(labels).toJSON();
        return new Basis({
            elements: elementsArrayWithIdsJSON as BasisSchema["elements"],
            coordinates: coordinatesArrayWithIdsJSON,
            units,
            cell,
            labels: labelsArrayWithIdsJSON,
        });
    }

    constructor(config: BasisConfig = Basis.defaultConfig) {
        super(config);
        const { elements, coordinates, units, labels } = config;
        this.elements = elements;
        this.labels = labels;
        this.coordinates = coordinates;
        this.units = units;
        this.cell = new Cell(config.cell);
        this.units = units || (ATOMIC_COORD_UNITS.crystal as BasisSchema["units"]);
        this._elements = Elements.fromObjects(elements);
        // @ts-ignore
        this._coordinates = Coordinates.fromObjects(coordinates);
        // @ts-ignore
        this._labels = ArrayWithIds.fromObjects(labels);
    }

    toJSON(exclude?: string[]): AnyObject {
        return {
            ...super.toJSON(exclude),
        };
    }

    get cellRounded() {
        return this.cell.vectorArraysRounded;
    }

    get elementsArray(): object[] {
        return this._elements.toJSON();
    }

    getElementsAsObject(): BasisSchema["elements"] {
        return this.elements;
    }

    get coordinatesAsArray() {
        return this._coordinates.values;
    }

    get isInCartesianUnits(): boolean {
        return this.units === ATOMIC_COORD_UNITS.cartesian;
    }

    get isInCrystalUnits(): boolean {
        return this.units === ATOMIC_COORD_UNITS.crystal;
    }

    toCartesian(): void {
        if (this.isInCartesianUnits) return;
        this._coordinates.mapArrayInPlace((point) => this.cell.convertPointToCartesian(point));
        this.units = ATOMIC_COORD_UNITS.cartesian as BasisSchema["units"];
    }

    toCrystal(): void {
        if (this.isInCrystalUnits) return;
        this._coordinates.mapArrayInPlace((point) => this.cell.convertPointToCrystal(point));
        this.units = ATOMIC_COORD_UNITS.crystal as BasisSchema["units"];
    }

    getElementByIndex(idx: number): string {
        return this._elements.getElementValueByIndex(idx) as string;
    }

    getCoordinateByIndex(idx: number): Coordinate {
        return Coordinate.fromArray(this._coordinates.getElementValueByIndex(idx) as PointSchema);
    }

    toStandardRepresentation() {
        this.toCrystal();
        this._coordinates.mapArrayInPlace((point) => point.map((x) => math.mod(x)) as PointSchema);
    }

    /** A representation where all coordinates are within 0 and 1 in crystal units */
    get standardRepresentation(): BasisSchema {
        const originalUnits = this.units;

        this.toStandardRepresentation();
        const result = this.toJSON();

        // preserve the original state
        if (originalUnits !== ATOMIC_COORD_UNITS.crystal) this.toCartesian();

        return result as unknown as BasisSchema;
    }

    /**
     * Add atom with a chemical element at coordinate.
     */
    // @ts-ignore
    addAtom({ element = "Si", coordinate = [0.5, 0.5, 0.5] }: Atom) {
        this._elements.addItem(element);
        // @ts-ignore
        this._coordinates.addItem(coordinate);
    }

    /**
     * Remove atom with a chemical element at coordinate either by passing the (element AND coordinate) OR id.
     */
    removeAtom({ element, coordinate, id }: Atom) {
        if (element && coordinate) {
            // @ts-ignore
            this._elements.removeItem(element, id);
            // @ts-ignore
            this._coordinates.removeItem(coordinate, id);
        }
    }

    /**
     * Unique names (symbols) of the chemical elements basis. E.g. `['Si', 'Li']`
     */
    get uniqueElements(): string[] {
        return _.unique(this._elements.values);
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
        return this._elements.values.map((element, idx) => {
            const coordinates = this.getCoordinateByIndex(idx);
            return [element, coordinates];
        });
    }

    /**
     * Returns a nested array with elements and their corresponding coordinates with labels
     * @example Output: [ ["Si", [0,0,0], ['1']], ["Si", [0.5,0.5,0.5]] , ['2']]
     */
    get elementsAndCoordinatesAndLabelsArray(): [string, Coordinate, string][] {
        return this._elements.values.map((element, idx) => {
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
        const clsInstance = this.clone();
        clsInstance.toStandardRepresentation();
        const standardRep = clsInstance.elementsAndCoordinatesAndLabelsArray.map((entry) => {
            const element = entry[0];
            const coordinate = entry[1];
            const atomicLabel = entry[2];
            const toleratedCoordinates = coordinate.value.value.map((x) =>
                math.round(x, HASH_TOLERANCE),
            );
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
                // @ts-ignore
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
            return `${element}${atomicLabel} ${coordinate.value.value
                .map((x) => s.sprintf("%14.9f", x).trim())
                .join(" ")}`;
        });
    }

    /**
     * @summary Returns number of atoms in material
     */
    get nAtoms(): number {
        return this._elements.values.length;
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
        return !this.cell.vectorArrays
            .map((vector, idx) => {
                return math.vEqualWithTolerance(vector, anotherBasisClsInstance.cell.vectorArrays[idx]);
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
        if (this._elements.values.length === 1) {
            const elementSymbol = this._elements.getElementValueByIndex(0);
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
        if (this._elements.values.length >= 2) {
            for (let i = 0; i < this._elements.values.length; i++) {
                for (let j = i + 1; j < this._elements.values.length; j++) {
                    const distance = math.vDist(
                        this._coordinates.getElementValueByIndex(i) as PointSchema,
                        this._coordinates.getElementValueByIndex(j) as PointSchema,
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
        const transposedBasisCoordinates = math.transpose(this._coordinates.values);
        const centerOfCoordinatesVectors = [];
        for (let i = 0; i < 3; i++) {
            const center = // @ts-ignore
                transposedBasisCoordinates[i].reduce((a, b) => a + b) /
                this._elements.values.length;
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
