import { ArrayWithIds } from "../abstract/array_with_ids";

export class AtomicConstraints {
    static fromArray(array) {
        return new AtomicConstraints({ values: array });
    }

    /**
     * Create atomic constraints.
     * @param {Object} config
     * @param {ArrayWithIds|Array} config.values
     */
    constructor(config = {}) {
        this.name = "atomic_constraints";
        this.values = new ArrayWithIds(config.values || []);
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

    getByIndex(idx) {
        return this.values.getArrayElementByIndex(idx) || [];
    }

    /**
     * Get constraints for an atom with index as string.
     * @param {Number} idx - atom index.
     * @param {Function} mapFn (OPTIONAL) - a function to be applied to each constraint. By default 0 or 1 is returned.
     * @return {String[]}
     */
    getAsStringByIndex(idx, mapFn = (boolean) => (boolean ? 1 : 0)) {
        return this.getByIndex(idx).map(mapFn).join(" ");
    }
}
