"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.ConstrainedBasis = void 0;
const underscore_string_1 = __importDefault(require("underscore.string"));
const array_with_ids_1 = require("../abstract/array_with_ids");
const constraints_1 = require("../constraints/constraints");
const basis_1 = require("./basis");
/**
 * @summary Extension of the Basis class able to deal with atomic constraints.
 * @extends Basis
 */
class ConstrainedBasis extends basis_1.Basis {
    /**
     * Create a an array with ids.
     * @param {Object} config
     * @param {ArrayWithIds|Array} config.constraints - atomic constraints.
     */
    constructor(config) {
        super(config);
        this._constraints = new array_with_ids_1.ArrayWithIds(config.constraints); // `constraints` is an Array with ids
    }
    get constraints() {
        return this._constraints.array;
    }
    set constraints(newConstraints) {
        this._constraints = new array_with_ids_1.ArrayWithIds(newConstraints);
    }
    getConstraintAsArray() {
        return this._constraints;
    }
    get AtomicConstraints() {
        return constraints_1.AtomicConstraints.fromArray(this.constraints);
    }
    /**
     * Serialize class instance to JSON.
     * @example As below:
         {
            ...Basis.toJSON(),
            "constraints": [
                {
                    "id" : 0,
                    "value" : [
                        1,
                        1,
                        1
                    ]
                },
            ]
         }
     */
    toJSON() {
        return {
            ...super.toJSON(),
            constraints: this._constraints.toJSON(),
        };
    }
    getConstraintByIndex(idx) {
        return this._constraints.getArrayElementByIndex(idx) || [];
    }
    /**
     * Helper function returning a nested array with [element, coordinates, constraints] as elements
     */
    get elementsCoordinatesConstraintsArray() {
        return this._elements.array.map((element, idx) => {
            const coordinates = this.getCoordinateByIndex(idx);
            const constraints = this.getConstraintByIndex(idx);
            return [element, coordinates, constraints];
        });
    }
    /**
     * Returns an array with atomic positions (with constraints) per atom stored as strings.
     * E.g., ``` ['Si  0 0 0  0 1 0', 'Li  0.5 0.5 0.5  1 0 1']```
     */
    get atomicPositionsWithConstraints() {
        return this.elementsCoordinatesConstraintsArray.map((entry) => {
            const element = entry[0];
            const coordinate = entry[1];
            const constraint = entry[2];
            return (underscore_string_1.default.sprintf("%-4s", element) +
                coordinate.map((x) => underscore_string_1.default.sprintf("%14.9f", x).trim()).join(" ") +
                " " +
                constraint.map((x) => (x ? 1 : 0)).join(" "));
        });
    }
}
exports.ConstrainedBasis = ConstrainedBasis;
//# sourceMappingURL=constrained_basis.js.map