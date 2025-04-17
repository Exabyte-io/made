import { AnyObject } from "@mat3ra/esse/dist/js/esse/types";
import {
    ArrayOf3NumberElementsSchema,
    LatticeExplicitUnit,
    MaterialSchema,
} from "@mat3ra/esse/dist/js/types";

export type MaterialJSON = MaterialSchema & AnyObject;

export type Vector = ArrayOf3NumberElementsSchema;
export type VectorsAsArray = [Vector, Vector, Vector];
export type latticeVectorUnits = Required<LatticeExplicitUnit>["units"];
