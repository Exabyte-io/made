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
        let _id;
        let _value;
        // if already passing a ScalarWithId => preserve original
        if (_.isObject(valueOrObject) && !_.isArray(valueOrObject)) {
            // NOTE - Arrays are Objects too
            _id = valueOrObject.id;
            _value = valueOrObject.value;
        } else {
            _id = id;
            _value = valueOrObject;
        }
        this.value = _value;
        this.id = _id;
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
