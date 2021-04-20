import s from "underscore.string";

import { ArrayWithIds } from "../abstract/array_with_ids";
import { AtomicConstraints } from "../constraints/constraints";
import { Basis } from "./basis";

/**
 * @summary Extension of the Basis class able to deal with atomic constraints.
 * @extends Basis
 */
export class ConstrainedBasis extends Basis {
    /**
     * Create a an array with ids.
     * @param {Object} config
     * @param {ArrayWithIds|Array} config.constraints - atomic constraints.
     */
    constructor(config) {
        super(config);
        this._constraints = new ArrayWithIds(config.constraints); // `constraints` is an Array with ids
    }

    get constraints() {
        return this._constraints;
    }

    set constraints(newConstraints) {
        this._constraints = new ArrayWithIds(newConstraints);
    }

    get AtomicConstraints() {
        return AtomicConstraints.fromArray(this.constraints.array);
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
            constraints: this.constraints.toJSON(),
        };
    }

    getConstraintByIndex(idx) {
        return this._constraints.getArrayElementByIndex(idx) || [];
    }

    /**
     * Helper function returning a nested array with [element, coordinates, constraints] as elements
     * @return {Array[]}
     */
    get elementsCoordinatesConstraintsArray() {
        const clsInstance = this;
        return this._elements.array.map((element, idx) => {
            const coordinates = clsInstance.getCoordinateByIndex(idx);
            const constraints = clsInstance.getConstraintByIndex(idx);
            return [element, coordinates, constraints];
        });
    }

    /**
     * Returns an array with atomic positions (with constraints) per atom stored as strings.
     * E.g., ``` ['Si  0 0 0  0 1 0', 'Li  0.5 0.5 0.5  1 0 1']```
     * @return {String[]}
     */
    get atomicPositionsWithConstraints() {
        return this.elementsCoordinatesConstraintsArray.map((entry) => {
            const element = entry[0];
            const coordinate = entry[1];
            const constraint = entry[2];
            return `${
                s.sprintf("%-4s", element) +
                coordinate.map((x) => s.sprintf("%14.9f", x).trim()).join(" ")
            } ${constraint.map((x) => (x ? 1 : 0)).join(" ")}`;
        });
    }
}
