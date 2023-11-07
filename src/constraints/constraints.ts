import { ArrayWithIds } from "../abstract/array_with_ids";
import { ObjectWithIdAndValue } from "../abstract/scalar_with_id";

export interface ConstraintValue extends Array<boolean> {
    0: boolean;
    1: boolean;
    2: boolean;
}

export type Constraint = ObjectWithIdAndValue<ConstraintValue>;

export class AtomicConstraints {
    name: string;

    values: ArrayWithIds<ConstraintValue>;

    static fromArray(array: ConstraintValue[]) {
        return new AtomicConstraints({ values: array });
    }

    /**
     * Create atomic constraints.
     * @param {Object} config
     * @param {ArrayWithIds|Array} config.values
     */
    constructor({ values }: { values?: ConstraintValue[] }) {
        this.name = "atomic_constraints";
        this.values = new ArrayWithIds(values || []);
    }

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
    toJSON() {
        return {
            name: this.name,
            values: this.values.toJSON(),
        };
    }

    getByIndex(idx: number): ConstraintValue {
        return this.values.getArrayElementByIndex(idx) || [];
    }

    /**
     * Get constraints for an atom with index as string.
     * @param idx - atom index.
     * @param mapFn (OPTIONAL) - a function to be applied to each constraint. By default 0 or 1 is returned.
     */
    getAsStringByIndex(idx: number, mapFn = (val: boolean): string => (val ? "1" : "0")): string {
        return this.getByIndex(idx).map(mapFn).join(" ");
    }
}
