import { LatticeImplicitSchema, LatticeTypeSchema } from "@mat3ra/esse/dist/js/types";
import { Vector, VectorsAsArray } from "./types";
export type Units = Required<LatticeImplicitSchema>["units"];
export interface FromVectorsProps {
    a: Vector;
    b: Vector;
    c: Vector;
    alat?: number;
    units?: Units["length"];
    type?: LatticeTypeSchema;
    skipRounding?: boolean;
}
export declare class LatticeBravais implements LatticeImplicitSchema {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
    units: Units;
    type: LatticeImplicitSchema["type"];
    /**
     * Create a Bravais lattice.
     */
    constructor(config: Partial<LatticeImplicitSchema>);
    static _roundValue(x: number): number;
    /**
     * Create a Bravais lattice from vectors.
     */
    static fromVectors({ a, b, c, alat, units, type, skipRounding, }: FromVectorsProps): LatticeBravais;
    /**
     * See fromVectors above.
     */
    static fromVectorArrays(array: VectorsAsArray, type: LatticeTypeSchema, skipRounding?: boolean): LatticeBravais;
    /**
     * Get the list of editable keys (eg. 'a', 'alpha') for the current lattice.
     * @return {Object}
     * @example {a: true, b: false, c: false, alpha: true, beta: false, gamma: false}
     */
    get editables(): {};
    /**
     * Serialize class instance to JSON.
     * @example As below:
         {
            "a" : 3.867,
            "b" : 3.867,
            "c" : 3.867,
            "alpha" : 60,
            "beta" : 60,
            "gamma" : 60,
            "units" : {
                "length" : "angstrom",
                "angle" : "degree"
            },
            "type" : "FCC"
         }
     */
    toJSON(): LatticeImplicitSchema;
}
