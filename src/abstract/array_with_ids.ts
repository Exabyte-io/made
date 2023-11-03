import _ from "underscore";

import { ScalarWithId } from "./scalar_with_id";

interface ObjectWithId {
    id: number;
    value?: object;
}

type Predicate = (o: object) => boolean;

type MapFunction = (value: object, index: number, array: object[]) => object;

/**
 * Helper class representing an ArrayWithIds. Used to explicitly track values assigned to atoms, for example.
 */
export class ArrayWithIds {
    array: object[];

    /**
     * Create a an array with ids.
     * @param {Array} array - Either regular array or ArrayWithIds (see @example above)
     */
    constructor(array: ObjectWithId[] = []) {
        if (!_.isArray(array))
            throw new Error("ArrayWithIds.constructor: pass array on initialization");
        // if passed an array with ids as config, only store the values in array
        this.array = array.sort((a, b) => a.id - b.id).map((element) => element.value || element);
    }

    /**
     * Serialize class instance to JSON.
     * @example [{"id" : 0, "value" : "Si" }, {"id" : 1, "value" : "Si" }]
     */
    toJSON() {
        // from ["a", "b"] to [{id: 0, value: "a"}, {id: 1, value: "b"}]
        return this.array.map((el, idx) => new ScalarWithId(el, idx).toJSON());
    }

    /**
     * Apply function fn to each element of the array and replace `array` with the result.
     * @param fn - The function to be applied to each array element.
     */
    mapArrayInPlace(fn: MapFunction) {
        if (!_.isFunction(fn))
            throw new Error("ArrayWithIds.mapArray: must pass function as argument");
        this.array = this.array.map((...args) => fn(...args));
    }

    getArrayElementByIndex(idx: number) {
        return this.array[idx];
    }

    /**
     * Get the index of the array element that passes the predicate.
     * @param {Function} predicate - The function to be applied to each array element.
     */
    getArrayIndexByPredicate(predicate: Predicate) {
        return this.array.findIndex((el) => predicate(el));
    }

    /**
     * Add an entity to array.
     * @param el - The entity to be added to array. If Object with 'value' key, its value will be added.
     */
    addElement(el: ObjectWithId) {
        if (el) this.array.push(el.value || el);
    }

    /**
     * Remove an entity to array. Either by passing the entity, or the corresponding index.
     * @param el - The entity to be added to array. If Object with 'value' key, its value will be added.
     * @param idx - The entity to be added to array. If Object with 'value' key, its value will be added.
     */
    removeElement(el: ObjectWithId, idx: number) {
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
