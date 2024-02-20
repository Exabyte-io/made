import { LatticeTypeSchema } from "@mat3ra/esse/lib/js/types";
export declare const PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS: {
    CUB: number[][];
    FCC: number[][];
    BCC: number[][];
    TET: number[][];
    BCT: number[][];
    ORC: number[][];
    ORCF: number[][];
    ORCI: number[][];
    ORCC: number[][];
    HEX: number[][];
    RHL: number[][];
    MCL: number[][];
    MCLC: number[][];
    TRI: number[][];
    TRIalt: number[][];
};
export declare const PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES: {
    [key in LatticeTypeSchema | "TRIalt"]: LatticeTypeSchema;
};
export declare function isConventionalCellSameAsPrimitiveForLatticeType(latticeType: LatticeTypeSchema): boolean;
