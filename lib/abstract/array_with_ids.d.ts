import { ObjectWithIdAndValue, ValueOrObject } from "./scalar_with_id";
type Predicate<T> = (o: T) => boolean;
type MapFunction<T> = (value: T, index: number, array: T[]) => T;
export declare function isArrayOfObjectsWithIdAndValue<T>(valueOrObjects: ValueOrObject<T>[]): valueOrObjects is ObjectWithIdAndValue<T>[];
/**
 * Helper class representing an ArrayWithIds. Used to explicitly track values assigned to atoms, for example.
 */
export declare class ArrayWithIds<T> {
    array: T[];
    /**
     * Create a an array with ids.
     * @param {Array} array - Either regular array or ArrayWithIds (see @example above)
     */
    constructor(array?: ObjectWithIdAndValue<T>[] | T[]);
    /**
     * Serialize class instance to JSON.
     * @example [{"id" : 0, "value" : "Si" }, {"id" : 1, "value" : "Si" }]
     */
    toJSON(): ObjectWithIdAndValue<T>[];
    /**
     * Apply function fn to each element of the array and replace `array` with the result.
     * @param fn - The function to be applied to each array element.
     */
    mapArrayInPlace(fn: MapFunction<T>): void;
    getArrayElementByIndex(idx: number): T;
    /**
     * Get the index of the array element that passes the predicate.
     * @param {Function} predicate - The function to be applied to each array element.
     */
    getArrayIndexByPredicate(predicate: Predicate<T>): number;
    /**
     * Add an entity to array.
     * @param el - The entity to be added to array. If Object with 'value' key, its value will be added.
     */
    addElement(el: ValueOrObject<T>): void;
    /**
     * Remove an entity to array. Either by passing the entity, or the corresponding index.
     * @param el - The entity to be added to array. If Object with 'value' key, its value will be added.
     * @param idx - The entity to be added to array. If Object with 'value' key, its value will be added.
     */
    removeElement(el: ValueOrObject<T> | null, idx?: number): void;
}
export {};
