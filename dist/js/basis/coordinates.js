"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Coordinates = exports.Coordinate = void 0;
const code_1 = require("@mat3ra/code");
const lodash_1 = require("lodash");
class Coordinate extends code_1.RoundedValueWithId {
    constructor({ value, id }) {
        super(id, value);
        this.value = value;
    }
    getValueAlongAxis(axis = "z") {
        const index = { x: 0, y: 1, z: 2 }[axis];
        return this.value[index];
    }
    prettyPrint(decimalPlaces = 9, padding = 14) {
        return this.value.map((x) => (0, lodash_1.padStart)(x.toFixed(decimalPlaces), padding)).join(" ");
    }
}
exports.Coordinate = Coordinate;
class Coordinates extends code_1.RoundedArrayWithIds {
    getValuesAlongAxis(axis = "z") {
        return this.values.map((coord) => {
            const coordinate = Coordinate.fromValueAndId(coord);
            return coordinate.getValueAlongAxis(axis);
        });
    }
    getMaxValueAlongAxis(axis = "z") {
        return Math.max(...this.getValuesAlongAxis(axis));
    }
    getMinValueAlongAxis(axis = "z") {
        return Math.min(...this.getValuesAlongAxis(axis));
    }
    getExtremumValueAlongAxis(extremum = "max", axis = "z") {
        return extremum === "max"
            ? this.getMaxValueAlongAxis(axis)
            : this.getMinValueAlongAxis(axis);
    }
}
exports.Coordinates = Coordinates;
