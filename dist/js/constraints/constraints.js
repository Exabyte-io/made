"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.AtomicConstraints = void 0;
const code_1 = require("@mat3ra/code");
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
