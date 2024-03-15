"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.ArrayWithIds = exports.isArrayOfObjectsWithIdAndValue = void 0;
const underscore_1 = __importDefault(require("underscore"));
const scalar_with_id_1 = require("./scalar_with_id");
function isArrayOfObjectsWithIdAndValue(valueOrObjects) {
    return (0, scalar_with_id_1.isObjectWithIdAndValue)(valueOrObjects[0]);
}
exports.isArrayOfObjectsWithIdAndValue = isArrayOfObjectsWithIdAndValue;
/**
 * Helper class representing an ArrayWithIds. Used to explicitly track values assigned to atoms, for example.
 */
class ArrayWithIds {
    /**
     * Create a an array with ids.
     * @param {Array} array - Either regular array or ArrayWithIds (see @example above)
     */
    constructor(array = []) {
        if (!underscore_1.default.isArray(array)) {
            throw new Error("ArrayWithIds.constructor: pass array on initialization");
        }
        // if passed an array with ids as config, only store the values in array
        if (isArrayOfObjectsWithIdAndValue(array)) {
            this.array = array.sort((a, b) => a.id - b.id).map((element) => element.value);
        }
        else {
            this.array = [...array];
        }
    }
    /**
     * Serialize class instance to JSON.
     * @example [{"id" : 0, "value" : "Si" }, {"id" : 1, "value" : "Si" }]
     */
    toJSON() {
        // from ["a", "b"] to [{id: 0, value: "a"}, {id: 1, value: "b"}]
        return this.array.map((el, idx) => new scalar_with_id_1.ScalarWithId(el, idx).toJSON());
    }
    /**
     * Apply function fn to each element of the array and replace `array` with the result.
     * @param fn - The function to be applied to each array element.
     */
    mapArrayInPlace(fn) {
        if (!underscore_1.default.isFunction(fn)) {
            throw new Error("ArrayWithIds.mapArray: must pass function as argument");
        }
        this.array = this.array.map(fn);
    }
    getArrayElementByIndex(idx) {
        return this.array[idx];
    }
    /**
     * Get the index of the array element that passes the predicate.
     * @param {Function} predicate - The function to be applied to each array element.
     */
    getArrayIndexByPredicate(predicate) {
        return this.array.findIndex((el) => predicate(el));
    }
    /**
     * Add an entity to array.
     * @param el - The entity to be added to array. If Object with 'value' key, its value will be added.
     */
    addElement(el) {
        const value = (0, scalar_with_id_1.isObjectWithIdAndValue)(el) ? el.value : el;
        if (el)
            this.array.push(value);
    }
    /**
     * Remove an entity to array. Either by passing the entity, or the corresponding index.
     * @param el - The entity to be added to array. If Object with 'value' key, its value will be added.
     * @param idx - The entity to be added to array. If Object with 'value' key, its value will be added.
     */
    removeElement(el, idx) {
        let _idx;
        if (idx === undefined) {
            _idx = this.array.findIndex((elm) => elm === el);
        }
        else {
            _idx = idx;
        }
        if (_idx !== undefined) {
            this.array.splice(_idx, 1);
        }
    }
}
exports.ArrayWithIds = ArrayWithIds;
//# sourceMappingURL=array_with_ids.js.map