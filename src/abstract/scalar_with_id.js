import _ from "underscore";

export class ScalarWithId {
    constructor(value, id = 0) {
        // if already passing a ScalarWithId => preserve original
        if (_.isObject(value) && !_.isArray(value)) { // Arrays are Objects too
            id = value.id;
            value = value.value;
        }
        this.value = value;
        this.id = id;
    }

    toJSON() {
        return {
            id: this.id,
            value: this.value,
        }
    }
}