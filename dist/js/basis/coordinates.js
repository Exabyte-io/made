"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.Coordinates = exports.Coordinate = void 0;
const code_1 = require("@mat3ra/code");
const lodash_1 = require("lodash");
const math_1 = __importDefault(require("../math"));
class Coordinate extends code_1.RoundedValueWithId {
    constructor({ value, id }) {
        super(id, value);
        this.value = value;
    }
    getValueAlongAxis(axis = "z") {
        const index = { x: 0, y: 1, z: 2 }[axis];
        return this.value[index];
    }
    translateByVector(vector) {
        this.value = this.value.map((v, i) => v + vector[i]);
        return this;
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
    translateByVector(vector) {
        this.mapArrayInPlace((x) => x.map((v, i) => v + vector[i]));
    }
    getCenterPoint() {
        const transposed = math_1.default.transpose(this.values);
        const center = [0, 0, 0];
        for (let i = 0; i < 3; i++) {
            const axisCoords = transposed[i];
            const sum = axisCoords.reduce((a, b) => a + b, 0);
            center[i] = math_1.default.precise(sum / this.values.length, 4);
        }
        return center;
    }
}
exports.Coordinates = Coordinates;
