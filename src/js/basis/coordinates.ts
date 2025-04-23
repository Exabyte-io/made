import { RoundedArrayWithIds, RoundedValueWithId } from "@mat3ra/code";
import {
    AtomicCoordinateSchema,
    Coordinate3DSchema,
    Matrix3X3Schema,
    Vector3DSchema,
} from "@mat3ra/esse/dist/js/types";
import { padStart } from "lodash";

import math from "../math";

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

    getCenterPoint(): Vector3DSchema {
        const transposed = math.transpose(this.values) as Matrix3X3Schema;
        const center: Vector3DSchema = [0, 0, 0];

        for (let i = 0; i < 3; i++) {
            const axisCoords = transposed[i] as Coordinate3DSchema;
            const sum = axisCoords.reduce((a, b) => a + b, 0);
            center[i] = math.precise(sum / this.values.length, 4);
        }

        return center;
    }
}
