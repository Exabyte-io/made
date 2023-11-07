import _ from "underscore";

export interface ObjectWithIdAndValue<T> {
    id: number;
    value: T;
}

export type ValueOrObject<T> = ObjectWithIdAndValue<T> | T;

export type ValueOrObjectArray<T> = ObjectWithIdAndValue<T>[] | T[];

export function isObjectWithIdAndValue<T>(
    valueOrObject: ValueOrObject<T>,
): valueOrObject is ObjectWithIdAndValue<T> {
    return Boolean(_.isObject(valueOrObject) && !_.isArray(valueOrObject) && valueOrObject.value);
}

/**
 * Helper class representing a scalar with an associated id.
 */
export class ScalarWithId<T> {
    id: number;

    value: T;

    /**
     * Create a an array with ids.
     * @param valueOrObject - a ScalarWithID, or any other type.
     * @param id - numerical id (Integer).
     */
    constructor(valueOrObject: ValueOrObject<T>, id = 0) {
        // if already passing a ScalarWithId => preserve original
        if (isObjectWithIdAndValue(valueOrObject)) {
            // NOTE - Arrays are Objects too
            this.id = valueOrObject.id;
            this.value = valueOrObject.value;
        } else {
            this.id = id;
            this.value = valueOrObject;
        }
    }

    /**
     * Serialize class instance to JSON.
     * @example {"id" : 0, "value" : "Si" }
     */
    toJSON(): ObjectWithIdAndValue<T> {
        return {
            id: this.id,
            value: this.value,
        };
    }
}
