import _ from "underscore";

import {
    isObjectWithIdAndValue,
    ObjectWithIdAndValue,
    ScalarWithId,
    ValueOrObject,
} from "./scalar_with_id";

type Predicate<T> = (o: T) => boolean;

type MapFunction<T> = (value: T, index: number, array: T[]) => T;

export function isArrayOfObjectsWithIdAndValue<T>(
    valueOrObjects: ValueOrObject<T>[],
): valueOrObjects is ObjectWithIdAndValue<T>[] {
    return isObjectWithIdAndValue(valueOrObjects[0]);
}

/**
 * Helper class representing an ArrayWithIds. Used to explicitly track values assigned to atoms, for example.
 */
export class ArrayWithIds<T> {
    array: T[];

    /**
     * Create a an array with ids.
     * @param {Array} array - Either regular array or ArrayWithIds (see @example above)
     */
    constructor(array: ObjectWithIdAndValue<T>[] | T[] = []) {
        if (!_.isArray(array)) {
            throw new Error("ArrayWithIds.constructor: pass array on initialization");
        }
        // if passed an array with ids as config, only store the values in array
        if (isArrayOfObjectsWithIdAndValue<T>(array)) {
            this.array = array.sort((a, b) => a.id - b.id).map((element) => element.value);
        } else {
            this.array = [...array];
        }
    }

    /**
     * Serialize class instance to JSON.
     * @example [{"id" : 0, "value" : "Si" }, {"id" : 1, "value" : "Si" }]
     */
    toJSON(): ObjectWithIdAndValue<T>[] {
        // from ["a", "b"] to [{id: 0, value: "a"}, {id: 1, value: "b"}]
        return this.array.map((el, idx) => new ScalarWithId<T>(el, idx).toJSON());
    }

    /**
     * Apply function fn to each element of the array and replace `array` with the result.
     * @param fn - The function to be applied to each array element.
     */
    mapArrayInPlace(fn: MapFunction<T>) {
        if (!_.isFunction(fn)) {
            throw new Error("ArrayWithIds.mapArray: must pass function as argument");
        }

        this.array = this.array.map(fn);
    }

    getArrayElementByIndex(idx: number) {
        return this.array[idx];
    }

    /**
     * Get the index of the array element that passes the predicate.
     * @param {Function} predicate - The function to be applied to each array element.
     */
    getArrayIndexByPredicate(predicate: Predicate<T>) {
        return this.array.findIndex((el) => predicate(el));
    }

    /**
     * Add an entity to array.
     * @param el - The entity to be added to array. If Object with 'value' key, its value will be added.
     */
    addElement(el: ValueOrObject<T>) {
        const value = isObjectWithIdAndValue(el) ? el.value : el;
        if (el) this.array.push(value);
    }

    /**
     * Remove an entity to array. Either by passing the entity, or the corresponding index.
     * @param el - The entity to be added to array. If Object with 'value' key, its value will be added.
     * @param idx - The entity to be added to array. If Object with 'value' key, its value will be added.
     */
    removeElement(el: ValueOrObject<T> | null, idx?: number) {
        let _idx;
        if (idx === undefined) {
            _idx = this.array.findIndex((elm) => elm === el);
        } else {
            _idx = idx;
        }
        if (_idx !== undefined) {
            this.array.splice(_idx, 1);
        }
    }
}
