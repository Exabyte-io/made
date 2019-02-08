import _ from "underscore";

class ScalarWithId {
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

export class ArrayWithIds {
    constructor(array = []) {
        if (!_.isArray(array)) throw new Error("ArrayWithIds.constructor: must pass array type on initialization");
        // if passed an array with ids as config, only store the values in array
        this.array = array.sort((a, b) => a.id - b.id).map(element => element.value || element);
    }

    // ["a", "b"] => [{id: 0, value: "a"}, {id: 1, value: "b"}]
    toJSON() {
        return this.array.map((el, idx) => new ScalarWithId(el, idx).toJSON());
    }

    // apply function fn to each element of the array and save the result in 'array'
    mapArrayInPlace(fn) {
        if (!_.isFunction(fn)) throw new Error("ArrayWithIds.mapArray: must pass function as argument");
        this.array = this.array.map((...args) => fn(...args));
    }

    getArrayElementByIndex(idx) {
        return this.array[idx];
    }

    getArrayIndexByPredicate(predicate) {
        return this.array.findIndex(el => predicate(el));
    }

    addElement(el) {
        el && this.array.push(el.value || el);
    }

    removeElement(el, idx) {
        if (idx === undefined) {
            idx = this.array.findIndex(elm => elm === el);
        }
        (idx !== undefined) && this.array.splice(idx, 1);
    }
}
