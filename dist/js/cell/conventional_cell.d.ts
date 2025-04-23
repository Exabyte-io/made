import { LatticeSchema, LatticeTypeEnum, Matrix3X3Schema } from "@mat3ra/esse/dist/js/types";
export declare const PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS: Record<string, Matrix3X3Schema>;
export declare const PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES: {
    [key in LatticeTypeEnum]: LatticeTypeEnum;
};
export declare function isConventionalCellSameAsPrimitiveForLatticeType(latticeType: LatticeSchema["type"]): boolean;
