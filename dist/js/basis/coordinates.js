"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.Coordinates = exports.Coordinate = void 0;
const code_1 = require("@mat3ra/code");
const underscore_string_1 = __importDefault(require("underscore.string"));
class Coordinate extends code_1.RoundedValueWithId {
    constructor({ value, id }) {
        super(id, value);
        this.value = value;
    }
    getValueAlongAxis(axis = "z") {
        const index = { x: 0, y: 1, z: 2 }[axis];
        return this.value[index];
    }
    prettyPrint(format = "%14.9f") {
        return this.value.map((x) => underscore_string_1.default.sprintf(format, x).trim()).join(" ");
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
