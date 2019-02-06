import s from "underscore.string";

import {Basis} from "./basis";
import {ArrayWithIds} from "../abstract/array_with_ids";
import {AtomicConstraints} from "../constraints/constraints";

/**
 * @summary Extension of the Basis class able to deal with atomic constraints
 * @param config {Object} Same as for Basis, plus constraints object
 * @param config.constraints {Array} Array of constraints containing AtomicConstraints.values
 */
export class ConstrainedBasis extends Basis {

    constructor(config) {
        super(config);
        this._constraints = new ArrayWithIds(config.constraints); // `constraints` is an Array with ids
    }

    get constraints() {
        return this._constraints
    }

    get AtomicConstraints() {
        return AtomicConstraints.fromArray(this.constraints.array);
    }

    set constraints(newConstraints) {
        this._constraints = new ArrayWithIds(newConstraints);
    }

    toJSON() {
        return {
            ...super.toJSON(),
            constraints: this.constraints.toJSON(),
        }
    }

    getConstraintByIndex(idx) {
        return this._constraints.getArrayElementByIndex(idx) || [];
    }

    get elementsCoordinatesConstraintsArray() {
        const clsInstance = this;
        return this._elements.array.map((element, idx) => {
            const coordinates = clsInstance.getCoordinateByIndex(idx);
            const constraints = clsInstance.getConstraintByIndex(idx);
            return [element, coordinates, constraints];
        });
    }

    /**
     * @summary Returns atomic positions array with constraints.
     * E.g., ``` ['Si 0 0 0', 'Li 0.5 0.5 0.5']```
     * @return {String[]}
     */
    get atomicPositionsWithConstraints() {
        const clsInstance = this;
        return this.elementsCoordinatesConstraintsArray.map((entry) => {
            const element = entry[0];
            const coordinate = entry[1];
            const constraint = entry[2];
            return s.sprintf('%-4s', element) +
                coordinate.map(x => s.sprintf('%14.9f', x).trim()).join(' ') +
                ' ' + constraint.map(x => x ? 1 : 0).join(' ');
        });
    }

}
