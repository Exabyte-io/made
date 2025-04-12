"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.LatticeBravais = void 0;
const constants_1 = __importDefault(require("../constants"));
const math_1 = __importDefault(require("../math"));
const types_1 = require("./types");
/*
 * @summary: class that holds parameters of a Bravais Lattice: a, b, c, alpha, beta, gamma + corresponding units.
 * When stored as class variables units for lengths are always "angstrom"s, angle - "degree"s
 */
class LatticeBravais {
    /**
     * Create a Bravais lattice.
     */
    constructor(config) {
        const { a = 1, // default lattice is cubic with unity in edge sizes
        b = a, c = a, alpha = 90, beta = alpha, gamma = alpha, 
        // if we do not know what lattice type this is => set to TRI
        type = "TRI", units = {
            length: "angstrom",
            angle: "degree",
        }, } = config;
        const k = constants_1.default.units.bohr === units.length ? constants_1.default.coefficients.BOHR_TO_ANGSTROM : 1;
        this.a = a * k;
        this.b = b * k;
        this.c = c * k;
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;
        this.type = type;
        this.units = units;
    }
    static _roundValue(x) {
        return math_1.default.precise(math_1.default.roundToZero(x));
    }
    /**
     * Create a Bravais lattice from vectors.
     */
    static fromVectors({ a, b, c, alat = 1, units = "angstrom", type = "TRI", skipRounding = false, }) {
        const roundValue = skipRounding ? (x) => x : this._roundValue;
        const config = {
            // @ts-ignore
            a: roundValue(math_1.default.vlen(a) * alat),
            // @ts-ignore
            b: roundValue(math_1.default.vlen(b) * alat),
            // @ts-ignore
            c: roundValue(math_1.default.vlen(c) * alat),
            alpha: roundValue(math_1.default.angle(b, c, "deg")),
            beta: roundValue(math_1.default.angle(a, c, "deg")),
            gamma: roundValue(math_1.default.angle(a, b, "deg")),
            // initially we do not know what lattice type this is => set to TRI
            type,
            units: {
                length: units,
                angle: "degree",
            },
        };
        // @ts-ignore
        return new this.prototype.constructor(config);
    }
    /**
     * See fromVectors above.
     */
    static fromVectorArrays(array, type, skipRounding = true) {
        return this.fromVectors({
            a: array[0],
            b: array[1],
            c: array[2],
            type,
            skipRounding, // do not round the values to avoid loosing precision by default
        });
    }
    /**
     * Get the list of editable keys (eg. 'a', 'alpha') for the current lattice.
     * @return {Object}
     * @example {a: true, b: false, c: false, alpha: true, beta: false, gamma: false}
     */
    get editables() {
        var _a;
        const object = {};
        const editablesList = (_a = types_1.LATTICE_TYPE_CONFIGS.find((entry) => entry.code === this.type)) === null || _a === void 0 ? void 0 : _a.editables;
        // ["a", "gamma"] => {a: true, gamma: true}
        if (editablesList) {
            editablesList.forEach((element) => {
                Object.assign(object, {
                    [element]: true,
                });
            });
        }
        return object;
    }
    /**
     * Serialize class instance to JSON.
     * @example As below:
         {
            "a" : 3.867,
            "b" : 3.867,
            "c" : 3.867,
            "alpha" : 60,
            "beta" : 60,
            "gamma" : 60,
            "units" : {
                "length" : "angstrom",
                "angle" : "degree"
            },
            "type" : "FCC"
         }
     */
    toJSON() {
        return {
            ...this,
            units: {
                length: this.units.length,
                angle: this.units.angle,
            },
        };
    }
}
exports.LatticeBravais = LatticeBravais;
