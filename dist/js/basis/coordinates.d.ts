import { RoundedArrayWithIds, RoundedValueWithId } from "@mat3ra/code";
import { AtomicCoordinateSchema, Coordinate3DSchema, Vector3DSchema } from "@mat3ra/esse/dist/js/types";
export type AtomicCoordinateValue = AtomicCoordinateSchema["value"];
type AxisType = "x" | "y" | "z";
export declare class Coordinate extends RoundedValueWithId<AtomicCoordinateValue> {
    value: AtomicCoordinateValue;
    constructor({ value, id }: AtomicCoordinateSchema);
    getValueAlongAxis(axis?: AxisType): number;
    translateByVector(vector: Vector3DSchema): Coordinate;
    prettyPrint(decimalPlaces?: number, padding?: number): string;
}
export declare class Coordinates extends RoundedArrayWithIds<Coordinate3DSchema> {
    getValuesAlongAxis(axis?: AxisType): number[];
    getMaxValueAlongAxis(axis?: AxisType): number;
    getMinValueAlongAxis(axis?: AxisType): number;
    getExtremumValueAlongAxis(extremum?: "max" | "min", axis?: AxisType): number;
    translateByVector(vector: Vector3DSchema): void;
}
export {};
