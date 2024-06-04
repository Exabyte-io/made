import { LatticeImplicitSchema, LatticeSchema, LatticeTypeExtendedSchema } from "@mat3ra/esse/dist/js/types";
import { Cell } from "../cell/cell";
import { FromVectorsProps, LatticeBravais } from "./lattice_bravais";
import { BravaisConfigProps, LatticeVectors } from "./lattice_vectors";
import { UnitCell } from "./unit_cell";
/**
 * Scaling factor used to calculate the new lattice size for non-periodic systems.
 * The scaling factor ensures that a non-periodic structure will have have a lattice greater than the structures size.
 */
export declare const nonPeriodicLatticeScalingFactor = 2;
export declare class Lattice extends LatticeBravais implements LatticeSchema {
    vectors: LatticeVectors;
    /**
     * Create a Lattice class from a config object.
     * @param {Object} config - Config object. See LatticeVectors.fromBravais.
     */
    constructor(config?: Partial<LatticeImplicitSchema>);
    /**
     * Create a Lattice class from a list of vectors.
     * @param {Object} config - Config object. See LatticeBravais.fromVectors.
     */
    static fromVectors(config: FromVectorsProps): Lattice;
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
            "type" : "FCC",
            "vectors" : {
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
        }
     */
    toJSON(skipRounding?: boolean): LatticeSchema;
    clone(extraContext: BravaisConfigProps): Lattice;
    /**
     * Get lattice vectors as a nested array
     */
    get vectorArrays(): [import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema];
    get Cell(): Cell;
    /**
     * Get a short label for the type of the lattice, eg. "MCLC".
     */
    get typeLabel(): string;
    /**
     * Get a short label for the extended type of the lattice, eg. "MCLC-5".
     */
    get typeExtended(): LatticeTypeExtendedSchema;
    /**
     * Calculate the volume of the lattice cell.
     */
    get volume(): number;
    static getDefaultPrimitiveLatticeConfigByType(latticeConfig: LatticeImplicitSchema): {
        a: number;
        b: number;
        c: number;
        alpha: number;
        beta: number;
        gamma: number;
        units: {
            length?: "angstrom" | "bohr" | undefined;
            angle?: "degree" | "radian" | undefined;
        };
        type: "CUB" | "BCC" | "FCC" | "TET" | "MCL" | "ORC" | "ORCC" | "ORCF" | "ORCI" | "HEX" | "BCT" | "TRI" | "MCLC" | "RHL";
    } & {
        a: number;
        b: number;
        c: number;
        alpha: number;
        beta: number;
        gamma: number;
    };
    get unitCell(): UnitCell;
    /**
     * Returns a string further used for the calculation of an unique hash.
     * @param isScaled - Whether to scale the vectors by the length of the first vector initially.
     */
    getHashString(isScaled?: boolean): string;
}
