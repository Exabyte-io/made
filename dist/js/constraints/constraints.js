"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.AtomicConstraints = exports.Constraint = void 0;
const code_1 = require("@mat3ra/code");
const lodash_1 = require("lodash");
class Constraint extends code_1.ValueWithId {
    constructor({ value, id }) {
        super({ id, value });
        this.value = value;
    }
    getValueAsString() {
        return this.value.map((x) => (x ? 1 : 0)).join(" ");
    }
    prettyPrint(element, coordinate, constraint) {
        return ((0, lodash_1.padEnd)(element, 4) +
            coordinate.prettyPrint() +
            " " +
            constraint.map((x) => (x ? 1 : 0)).join(" "));
    }
}
exports.Constraint = Constraint;
class AtomicConstraints extends code_1.ArrayWithIds {
    /**
     * Get constraints for an atom with index as string.
     * @param idx - atom index.
     * @param mapFn (OPTIONAL) - a function to be applied to each constraint. By default 0 or 1 is returned.
     */
    getAsStringByIndex(idx, mapFn = (val) => (val ? "1" : "0")) {
        const constraints = this.getElementValueByIndex(idx);
        return constraints ? constraints.map(mapFn).join(" ") : "";
    }
}
exports.AtomicConstraints = AtomicConstraints;
