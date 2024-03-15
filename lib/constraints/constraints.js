"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.AtomicConstraints = void 0;
const array_with_ids_1 = require("../abstract/array_with_ids");
class AtomicConstraints {
    static fromArray(array) {
        return new AtomicConstraints({ values: array });
    }
    /**
     * Create atomic constraints.
     * @param {Object} config
     * @param {ArrayWithIds|Array} config.values
     */
    constructor({ values }) {
        this.name = "atomic_constraints";
        this.values = new array_with_ids_1.ArrayWithIds(values || []);
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
     * @param idx - atom index.
     * @param mapFn (OPTIONAL) - a function to be applied to each constraint. By default 0 or 1 is returned.
     */
    getAsStringByIndex(idx, mapFn = (val) => (val ? "1" : "0")) {
        return this.getByIndex(idx).map(mapFn).join(" ");
    }
}
exports.AtomicConstraints = AtomicConstraints;
//# sourceMappingURL=constraints.js.map