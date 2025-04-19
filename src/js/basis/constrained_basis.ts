// @ts-nocheck
import { ArrayWithIds, ValueWithId } from "@mat3ra/code";
import { AtomicConstraints as ConstraintsSchema, BasisSchema } from "@mat3ra/esse/dist/js/types";
import s from "underscore.string";

import { AtomicConstraints, Constraint, ConstraintValue } from "../constraints/constraints";
import { Basis, BasisConfig } from "./basis";
import { Coordinate } from "./types";

export interface ConstrainedBasisProps extends BasisConfig {
    constraints: Constraint[];
}

export interface ConstrainedBasisJSON extends BasisSchema {
    constraints: ValueWithId<ConstraintValue>[];
}

/**
 * @summary Extension of the Basis class able to deal with atomic constraints.
 * @extends Basis
 */
export class ConstrainedBasis extends Basis {
    _constraints: ArrayWithIds<ConstraintValue>;

    constraints: Constraint[];

    /**
     * Create a an array with ids.
     * @param {Object} config
     * @param {ArrayWithIds|Array} config.constraints - atomic constraints.
     */
    constructor(config: ConstrainedBasisProps) {
        super(config);
        this.constraints = config.constraints;
        this._constraints = new ArrayWithIds<ConstraintValue>(config.constraints); // `constraints` is an Array with ids
    }

    getConstraintAsArray() {
        return this._constraints;
    }

    get AtomicConstraints() {
        return AtomicConstraints.fromArray(this.constraints);
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
    toJSON(): ConstrainedBasisJSON {
        return {
            ...super.toJSON(),
            constraints: this._constraints.toJSON(),
        };
    }

    getConstraintByIndex(idx: number): ConstraintValue {
        return this._constraints.getElementValueByIndex(idx) || [false, false, false];
    }

    /**
     * Helper function returning a nested array with [element, coordinates, constraints] as elements
     */
    get elementsCoordinatesConstraintsArray(): [string, Coordinate, ConstraintValue, string][] {
        return this._elements.values.map((element: any, idx: number) => {
            const coordinates = this.getCoordinateByIndex(idx);
            const constraints = this.getConstraintByIndex(idx);
            const atomicLabel = this.atomicLabelsArray[idx];
            return [element, coordinates, constraints, atomicLabel];
        });
    }

    /**
     * Returns an array with atomic positions (with constraints) per atom stored as strings.
     * E.g., ``` ['Si  0 0 0  0 1 0', 'Li  0.5 0.5 0.5  1 0 1']```
     */
    get atomicPositionsWithConstraints(): string[] {
        return this.elementsCoordinatesConstraintsArray.map((entry) => {
            const element = entry[0] + entry[3]; // element with label, Fe1
            const coordinate = entry[1];
            const constraint = entry[2];
            return (
                s.sprintf("%-4s", element) +
                coordinate.map((x) => s.sprintf("%14.9f", x).trim()).join(" ") +
                " " +
                constraint.map((x) => (x ? 1 : 0)).join(" ")
            );
        });
    }
}
