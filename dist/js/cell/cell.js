"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.Cell = void 0;
const code_1 = require("@mat3ra/code");
const math_1 = require("@mat3ra/code/dist/js/math");
const constants_1 = __importDefault(require("../constants"));
const MATRIX = math_1.math.matrix;
const MULT = math_1.math.multiply;
const INV = math_1.math.inv;
// @ts-ignore
const MATRIX_MULT = (...args) => MULT(...args.map((x) => MATRIX(x))).toArray();
class Cell {
    constructor(config = Cell.defaultConfig) {
        this.alat = 1;
        this.units = "angstrom";
        const { a, b, c } = config;
        this.a = a;
        this.b = b;
        this.c = c;
    }
    static fromVectorsArray(vectors) {
        return new Cell({
            a: vectors[0],
            b: vectors[1],
            c: vectors[2],
        });
    }
    get _a() {
        return new Cell.RoundedVector3DClassReference(this.a);
    }
    get _b() {
        return new Cell.RoundedVector3DClassReference(this.b);
    }
    get _c() {
        return new Cell.RoundedVector3DClassReference(this.c);
    }
    get vectorArrays() {
        return [this._a.value, this._b.value, this._c.value];
    }
    get vectorArraysRounded() {
        return [this._a.value_rounded, this._b.value_rounded, this._c.value_rounded];
    }
    get volume() {
        return math_1.math.det(this.vectorArrays);
    }
    get volumeRounded() {
        return math_1.math.roundArrayOrNumber(this.volume, Cell.roundPrecision);
    }
    clone() {
        return new Cell({ a: this.a, b: this.b, c: this.c });
    }
    cloneAndScaleByMatrix(matrix) {
        const newCell = this.clone();
        newCell.scaleByMatrix(matrix);
        return newCell;
    }
    convertPointToCartesian(point) {
        return MULT(point, this.vectorArrays);
    }
    convertPointToCrystal(point) {
        return MULT(point, INV(this.vectorArrays));
    }
    isPointInsideCellCartesian(point) {
        const { tolerance } = this.constructor;
        return (this.convertPointToCrystal(point)
            .map((c) => math_1.math.isBetweenZeroInclusiveAndOne(c, tolerance))
            // @ts-ignore
            .reduce((a, b) => a && b));
    }
    // eslint-disable-next-line class-methods-use-this
    isPointInsideCellCrystal(point) {
        const { tolerance } = this.constructor;
        return (point
            .map((c) => math_1.math.isBetweenZeroInclusiveAndOne(c, tolerance))
            // @ts-ignore
            .reduce((a, b) => a && b));
    }
    isPointInsideCell(point, useCrystal = false) {
        if (useCrystal) {
            return this.isPointInsideCellCrystal(point);
        }
        return this.isPointInsideCellCartesian(point);
    }
    getMostCollinearVectorIndex(testVector) {
        const angles = this.vectorArrays.map((v) => math_1.math.angleUpTo90(v, testVector, "deg"));
        return angles.findIndex((el) => el === math_1.math.min(angles));
    }
    scaleByMatrix(matrix) {
        [this.a, this.b, this.c] = MATRIX_MULT(matrix, this.vectorArrays);
    }
    toJSON() {
        return {
            a: this._a.toJSON(),
            b: this._b.toJSON(),
            c: this._c.toJSON(),
            alat: this.alat,
            units: this.units,
        };
    }
}
exports.Cell = Cell;
Cell.RoundedVector3DClassReference = code_1.RoundedVector3D;
Cell.roundPrecision = 9;
Cell.tolerance = constants_1.default.tolerance.length; // in crystal coordinates
Cell.defaultConfig = {
    a: [1.0, 0, 0],
    b: [0, 1.0, 0],
    c: [0, 0, 1.0],
};
