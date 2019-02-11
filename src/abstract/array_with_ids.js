import _ from "underscore";

import {ScalarWithId} from "./scalar_with_id";

/**
 * A helper class representing an ArrayWithIds.
 * @example [{"id" : 0, "value" : "Si" }, {"id" : 1, "value" : "Si" }]
 */
export class ArrayWithIds {
    /**
     * Create a an array with ids.
     * @param {Array} array - Either regular array or ArrayWithIds (see @example above)
     */
    constructor(array = []) {
        if (!_.isArray(array)) throw new Error("ArrayWithIds.constructor: pass array on initialization");
        // if passed an array with ids as config, only store the values in array
        this.array = array.sort((a, b) => a.id - b.id).map(element => element.value || element);
    }

    toJSON() {
        // from ["a", "b"] to [{id: 0, value: "a"}, {id: 1, value: "b"}]
        return this.array.map((el, idx) => new ScalarWithId(el, idx).toJSON());
    }

    /**
     * Apply function fn to each element of the array and save the result in 'array'
     * @param {Function} fn - The function to be applied to each array element.
     */
    mapArrayInPlace(fn) {
        if (!_.isFunction(fn)) throw new Error("ArrayWithIds.mapArray: must pass function as argument");
        this.array = this.array.map((...args) => fn(...args));
    }

    getArrayElementByIndex(idx) {return this.array[idx]}

    /**
     * Get the index of the array element that passes the predicate.
     * @param {Function} predicate - The function to be applied to each array element.
     */
    getArrayIndexByPredicate(predicate) {return this.array.findIndex(el => predicate(el))}

    /**
     * Add an entity to array.
     * @param {Any} el - The entity to be added to array. If Object with 'value' key, its value will be added.
     */
    addElement(el) {el && this.array.push(el.value || el)}

    /**
     * Remove an entity to array. Either by passing the entity, or the corresponding index.
     * @param {Any} el - The entity to be added to array. If Object with 'value' key, its value will be added.
     * @param {Any} idx - The entity to be added to array. If Object with 'value' key, its value will be added.
     */
    removeElement(el, idx) {
        if (idx === undefined) idx = this.array.findIndex(elm => elm === el);
        (idx !== undefined) && this.array.splice(idx, 1);
    }
}
