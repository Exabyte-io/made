"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.ScalarWithId = exports.isObjectWithIdAndValue = void 0;
const underscore_1 = __importDefault(require("underscore"));
function isObjectWithIdAndValue(valueOrObject) {
    return Boolean(underscore_1.default.isObject(valueOrObject) && !underscore_1.default.isArray(valueOrObject) && valueOrObject.value);
}
exports.isObjectWithIdAndValue = isObjectWithIdAndValue;
/**
 * Helper class representing a scalar with an associated id.
 */
class ScalarWithId {
    /**
     * Create a an array with ids.
     * @param valueOrObject - a ScalarWithID, or any other type.
     * @param id - numerical id (Integer).
     */
    constructor(valueOrObject, id = 0) {
        // if already passing a ScalarWithId => preserve original
        if (isObjectWithIdAndValue(valueOrObject)) {
            // NOTE - Arrays are Objects too
            this.id = valueOrObject.id;
            this.value = valueOrObject.value;
        }
        else {
            this.id = id;
            this.value = valueOrObject;
        }
    }
    /**
     * Serialize class instance to JSON.
     * @example {"id" : 0, "value" : "Si" }
     */
    toJSON() {
        return {
            id: this.id,
            value: this.value,
        };
    }
}
exports.ScalarWithId = ScalarWithId;
