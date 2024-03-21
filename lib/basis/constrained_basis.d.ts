import { ArrayWithIds } from "../abstract/array_with_ids";
import { ObjectWithIdAndValue } from "../abstract/scalar_with_id";
import { AtomicConstraints, Constraint, ConstraintValue } from "../constraints/constraints";
import { Basis, BasisProps, BasisSchema } from "./basis";
import { Coordinate } from "./types";
export interface ConstrainedBasisProps extends BasisProps {
    constraints: Constraint[];
}
export interface ConstrainedBasisJSON extends BasisSchema {
    constraints: ObjectWithIdAndValue<ConstraintValue>[];
}
/**
 * @summary Extension of the Basis class able to deal with atomic constraints.
 * @extends Basis
 */
export declare class ConstrainedBasis extends Basis {
    _constraints: ArrayWithIds<ConstraintValue>;
    /**
     * Create a an array with ids.
     * @param {Object} config
     * @param {ArrayWithIds|Array} config.constraints - atomic constraints.
     */
    constructor(config: ConstrainedBasisProps);
    get constraints(): ConstraintValue[];
    set constraints(newConstraints: ConstraintValue[]);
    getConstraintAsArray(): ArrayWithIds<ConstraintValue>;
    get AtomicConstraints(): AtomicConstraints;
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
    toJSON(): ConstrainedBasisJSON;
    getConstraintByIndex(idx: number): ConstraintValue;
    /**
     * Helper function returning a nested array with [element, coordinates, constraints] as elements
     */
    get elementsCoordinatesConstraintsArray(): [string, Coordinate, ConstraintValue, string][];
    /**
     * Returns an array with atomic positions (with constraints) per atom stored as strings.
     * E.g., ``` ['Si  0 0 0  0 1 0', 'Li  0.5 0.5 0.5  1 0 1']```
     */
    get atomicPositionsWithConstraints(): string[];
}
