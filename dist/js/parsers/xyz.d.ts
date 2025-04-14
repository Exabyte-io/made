import { MaterialSchema } from "@mat3ra/esse/dist/js/types";
import { ConstrainedBasis } from "../basis/constrained_basis";
import { Constraint, ConstraintValue } from "../constraints/constraints";
import { Vector } from "../lattice/types";
import { CombinatorialBasis } from "./xyz_combinatorial_basis";
/**
 * Validates that passed string is well-formed XYZ file.
 */
export declare function validate(xyzTxt: string): void;
export interface ParsedObject {
    element: string;
    coordinates: Vector;
    constraints: ConstraintValue;
    label?: number;
}
export interface BasisConfig {
    elements: {
        id: number;
        value: string;
    }[];
    labels?: {
        id: number;
        value: number;
    }[];
    coordinates: {
        id: number;
        value: Vector;
    }[];
    units: string;
    cell: Vector[];
    constraints: Constraint[];
}
/**
 * Parse XYZ text for basis.
 * @param txt Text
 * @param units Coordinate units
 * @param cell Basis Cell
 */
declare function toBasisConfig(txt: string, units?: string, cell?: [import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/dist/js/types").ArrayOf3NumberElementsSchema]): BasisConfig;
/**
 * Create XYZ from Basis class instance.
 * @param basisClsInstance Basis class instance.
 * @param printFormat Output format for coordinates.
 * @param skipRounding Whether to round the numbers (ie. to avoid negative zeros).
 * @return Basis string in XYZ format
 */
declare function fromBasis(basisClsInstance: ConstrainedBasis, printFormat?: string, skipRounding?: boolean): string;
/**
 * Create XYZ from Material class instance (or its JSON config).
 * @param materialOrConfig Material.
 * @param fractional Coordinate units as fractional.
 * @return Class Instance
 */
declare function fromMaterial(materialOrConfig: MaterialSchema, fractional?: boolean): string;
declare const _default: {
    validate: typeof validate;
    fromMaterial: typeof fromMaterial;
    toBasisConfig: typeof toBasisConfig;
    fromBasis: typeof fromBasis;
    CombinatorialBasis: typeof CombinatorialBasis;
};
export default _default;
