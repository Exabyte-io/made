import { LatticeExplicitUnit, LatticeImplicitSchema } from "@mat3ra/esse/dist/js/types";
import { Vector } from "./types";
type RequiredLatticeExplicitUnit = Required<LatticeExplicitUnit>;
export interface BravaisConfigProps extends Partial<LatticeImplicitSchema> {
    isConventional?: boolean;
}
export declare class LatticeVectors implements RequiredLatticeExplicitUnit {
    a: RequiredLatticeExplicitUnit["a"];
    b: RequiredLatticeExplicitUnit["b"];
    c: RequiredLatticeExplicitUnit["c"];
    alat: RequiredLatticeExplicitUnit["alat"];
    units: RequiredLatticeExplicitUnit["units"];
    /**
     * Create a Bravais lattice.
     */
    constructor(config: LatticeExplicitUnit);
    static _roundValue(arr: number[]): number[];
    static fromBravais({ a, // default lattice is cubic with unity in edge sizes
    b, c, alpha, beta, gamma, units, type, isConventional, }: BravaisConfigProps): LatticeVectors;
    get vectorArrays(): [Vector, Vector, Vector];
    /**
     * Serialize class instance to JSON.
     * @example As below:
     {
            "a" : [
                3.34892,
                0,
                1.9335
            ],
            "b" : [
                1.116307,
                3.157392,
                1.9335
            ],
            "c" : [
                0,
                0,
                3.867
            ],
            "alat" : 1,
            "units" : "angstrom"
        }
     */
    toJSON(): RequiredLatticeExplicitUnit;
}
export {};
