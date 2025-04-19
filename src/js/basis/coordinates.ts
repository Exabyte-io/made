import { RoundedArrayWithIds, RoundedValueWithId, RoundedVector3D } from "@mat3ra/code";
import { PointSchema } from "@mat3ra/esse/dist/js/types";

export class Coordinate extends RoundedValueWithId<RoundedVector3D> {
    value: RoundedVector3D;

    constructor(id: number, value: RoundedVector3D) {
        super(id, value);
        this.value = value;
    }

    static fromArray(value: number[], id = 0): Coordinate {
        return new Coordinate(id, new RoundedVector3D(value));
    }

    getValueAlongAxis(axis: "x" | "y" | "z" = "z"): number {
        return this.value[axis];
    }
}

export class Coordinates extends RoundedArrayWithIds<PointSchema> {
    getValuesAlongAxis(axis: "x" | "y" | "z" = "z"): number[] {
        return this.values.map((coord) => Coordinate.fromArray(coord).getValueAlongAxis(axis));
    }

    getMaxValueAlongAxis(axis: "x" | "y" | "z" = "z"): number {
        return Math.max(...this.getValuesAlongAxis(axis));
    }

    getMinValueAlongAxis(axis: "x" | "y" | "z" = "z"): number {
        return Math.min(...this.getValuesAlongAxis(axis));
    }

    getExtremumValueAlongAxis(
        extremum: "max" | "min" = "max",
        axis: "x" | "y" | "z" = "z",
    ): number {
        return extremum === "max"
            ? this.getMaxValueAlongAxis(axis)
            : this.getMinValueAlongAxis(axis);
    }
}
