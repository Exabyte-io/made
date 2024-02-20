import { ArrayWithIds } from "../abstract/array_with_ids";
import { ObjectWithIdAndValue } from "../abstract/scalar_with_id";
export interface ConstraintValue extends Array<boolean> {
    0: boolean;
    1: boolean;
    2: boolean;
}
export type Constraint = ObjectWithIdAndValue<ConstraintValue>;
export declare class AtomicConstraints {
    name: string;
    values: ArrayWithIds<ConstraintValue>;
    static fromArray(array: ConstraintValue[]): AtomicConstraints;
    /**
     * Create atomic constraints.
     * @param {Object} config
     * @param {ArrayWithIds|Array} config.values
     */
    constructor({ values }: {
        values?: ConstraintValue[];
    });
    /**
     * @example As below:
        [
            {
                "id" : 0,
                "value" : [
                    1,
                    1,
                    1
                ]
            }
        ]
     */
    toJSON(): {
        name: string;
        values: ObjectWithIdAndValue<ConstraintValue>[];
    };
    getByIndex(idx: number): ConstraintValue;
    /**
     * Get constraints for an atom with index as string.
     * @param idx - atom index.
     * @param mapFn (OPTIONAL) - a function to be applied to each constraint. By default 0 or 1 is returned.
     */
    getAsStringByIndex(idx: number, mapFn?: (val: boolean) => string): string;
}
