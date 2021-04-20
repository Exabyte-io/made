import _ from "underscore";

/**
 * Helper class representing a scalar with an associated id.
 */
export class ScalarWithId {
    /**
     * Create a an array with ids.
     * @param {Any} valueOrObject - a ScalarWithID, or any other type.
     * @param {Number} id - numerical id (Integer).
     */
    constructor(valueOrObject, id = 0) {
        // if already passing a ScalarWithId => preserve original
        if (_.isObject(valueOrObject) && !_.isArray(valueOrObject)) {
            // NOTE - Arrays are Objects too
            id = valueOrObject.id;
            valueOrObject = valueOrObject.value;
        }
        this.value = valueOrObject;
        this.id = id;
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
