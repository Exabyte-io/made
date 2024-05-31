import MathModule from "@mat3ra/code/dist/js/math";
import { Coordinate } from "../basis/types";
import { Vector, VectorsAsArray } from "../lattice/types";
type Point = Coordinate | MathModule.Matrix | MathModule.MathType;
export declare class Cell {
    tolerance: number;
    vector1: Vector;
    vector2: Vector;
    vector3: Vector;
    /**
     * Create a cell.
     * @param nestedArray {Number[][]} is an array of cell vectors in cartesian Angstrom units.
     */
    constructor(nestedArray: VectorsAsArray);
    /**
     * Get cell vectors as (a nested) array.
     * @example [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
     */
    get vectorsAsArray(): VectorsAsArray;
    clone(): Cell;
    cloneAndScaleByMatrix(matrix: number[][]): Cell;
    /**
     * Convert a point (in crystal coordinates) to cartesian.
     */
    convertPointToCartesian(point: Point): MathModule.MathType;
    /**
     * Convert a point (in cartesian coordinates) to crystal (fractional).
     */
    convertPointToFractional(point: Point): Coordinate;
    /**
     * Check whether a point is inside the cell.
     * @param point - the point to conduct the check for.
     * @param tolerance - numerical tolerance.
     */
    isPointInsideCell(point: Point, tolerance?: number): boolean;
    /**
     * Returns the index of the cell vector, most collinear with the testVector.
     * @param testVector
     */
    getMostCollinearVectorIndex(testVector: Vector): number;
    /**
     * Scale this cell by right-multiplying it to a matrix (nested array)
     */
    scaleByMatrix(matrix: number[][]): void;
}
export {};
