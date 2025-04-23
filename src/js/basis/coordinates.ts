import { RoundedArrayWithIds, RoundedValueWithId } from "@mat3ra/code";
import {
    AtomicCoordinateSchema,
    Coordinate3DSchema,
    Vector3DSchema,
} from "@mat3ra/esse/dist/js/types";
import { padStart } from "lodash";

export type AtomicCoordinateValue = AtomicCoordinateSchema["value"];

type AxisType = "x" | "y" | "z";

export class Coordinate extends RoundedValueWithId<AtomicCoordinateValue> {
    value: AtomicCoordinateValue;

    constructor({ value, id }: AtomicCoordinateSchema) {
        super(id, value);
        this.value = value;
    }

    getValueAlongAxis(axis: AxisType = "z"): number {
        const index = { x: 0, y: 1, z: 2 }[axis];
        return this.value[index] as number;
    }

    translateByVector(vector: Vector3DSchema): Coordinate {
        this.value = this.value.map((v, i) => v + vector[i]) as AtomicCoordinateValue;
        return this;
    }

    prettyPrint(decimalPlaces = 9, padding = 14): string {
        return this.value.map((x: number) => padStart(x.toFixed(decimalPlaces), padding)).join(" ");
    }
}

export class Coordinates extends RoundedArrayWithIds<Coordinate3DSchema> {
    getValuesAlongAxis(axis: AxisType = "z"): number[] {
        return this.values.map((coord) => {
            const coordinate = Coordinate.fromValueAndId(coord);
            return coordinate.getValueAlongAxis(axis);
        });
    }

    getMaxValueAlongAxis(axis: AxisType = "z"): number {
        return Math.max(...this.getValuesAlongAxis(axis));
    }

    getMinValueAlongAxis(axis: AxisType = "z"): number {
        return Math.min(...this.getValuesAlongAxis(axis));
    }

    getExtremumValueAlongAxis(extremum: "max" | "min" = "max", axis: AxisType = "z"): number {
        return extremum === "max"
            ? this.getMaxValueAlongAxis(axis)
            : this.getMinValueAlongAxis(axis);
    }

    translateByVector(vector: Vector3DSchema): void {
        this.mapArrayInPlace((x) => x.map((v, i) => v + vector[i]) as Coordinate3DSchema);
    }
}
