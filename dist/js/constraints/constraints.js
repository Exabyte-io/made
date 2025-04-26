"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.AtomicConstraints = exports.Constraint = void 0;
const code_1 = require("@mat3ra/code");
class Constraint extends code_1.ValueWithId {
    constructor({ value, id }) {
        super({ id, value });
        this.value = value;
    }
    getValueAsString() {
        return this.value.map((x) => (x ? 1 : 0)).join(" ");
    }
    prettyPrint() {
        return this.value.map((x) => (x ? 1 : 0)).join(" ");
    }
    // By default, the constraint is unconstrained if all values are 1/true.
    isUnconstrained() {
        return this.value.every((val) => val);
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
    get areUnconstrained() {
        return this.values.every((constraint) => {
            const _constraint = Constraint.fromValueAndId(constraint);
            return _constraint.isUnconstrained();
        });
    }
}
exports.AtomicConstraints = AtomicConstraints;
