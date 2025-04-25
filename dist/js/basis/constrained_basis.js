"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.ConstrainedBasis = void 0;
const constraints_1 = require("../constraints/constraints");
const basis_1 = require("./basis");
/**
 * @summary Extension of the Basis class able to deal with atomic constraints.
 * @extends Basis
 */
class ConstrainedBasis extends basis_1.Basis {
    constructor(config) {
        super(config);
        const { constraints } = config;
        this._constraints = constraints_1.AtomicConstraints.fromObjects(constraints || []); // `constraints` is an Array with ids
    }
    static fromElementsCoordinatesAndConstraints(config) {
        const basisConfig = this._convertValuesToConfig(config);
        const constraints = constraints_1.AtomicConstraints.fromValues(config.constraints);
        return new this({
            ...basisConfig,
            constraints: constraints.toJSON(),
        });
    }
    get constraints() {
        return this._constraints.toJSON();
    }
    set constraints(constraints) {
        this._constraints = constraints_1.AtomicConstraints.fromObjects(constraints);
    }
    get AtomicConstraints() {
        return constraints_1.AtomicConstraints.fromObjects(this.constraints);
    }
    toJSON() {
        return {
            ...super.toJSON(),
            constraints: this.constraints,
        };
    }
    getConstraintByIndex(idx) {
        return this._constraints.getElementValueByIndex(idx) || [false, false, false];
    }
    /**
     * Helper function returning a nested array with [element, coordinates, constraints] as elements
     */
    get elementsCoordinatesConstraintsArray() {
        return this._elements.values.map((element, idx) => {
            const coordinate = this.getCoordinateByIndex(idx);
            const constraint = this.getConstraintByIndex(idx);
            const atomicLabel = this.atomicLabelsArray[idx];
            return [element, coordinate, constraint, atomicLabel];
        });
    }
    /**
     * Returns an array with atomic positions (with constraints) per atom stored as strings.
     * E.g., ``` ['Si  0 0 0  0 1 0', 'Li  0.5 0.5 0.5  1 0 1']```
     */
    get atomicPositionsWithConstraints() {
        return this.elementsCoordinatesConstraintsArray.map(([element, coordinate, constraint, label]) => {
            const fullElement = element + label; // e.g., Fe1
            return new constraints_1.Constraint({ id: 0, value: constraint }).prettyPrint(fullElement, coordinate, constraint);
        });
    }
}
exports.ConstrainedBasis = ConstrainedBasis;
