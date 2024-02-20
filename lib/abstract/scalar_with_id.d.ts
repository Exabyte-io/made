export interface ObjectWithIdAndValue<T> {
    id: number;
    value: T;
}
export type ValueOrObject<T> = ObjectWithIdAndValue<T> | T;
export type ValueOrObjectArray<T> = ObjectWithIdAndValue<T>[] | T[];
export declare function isObjectWithIdAndValue<T>(valueOrObject: ValueOrObject<T>): valueOrObject is ObjectWithIdAndValue<T>;
/**
 * Helper class representing a scalar with an associated id.
 */
export declare class ScalarWithId<T> {
    id: number;
    value: T;
    /**
     * Create a an array with ids.
     * @param valueOrObject - a ScalarWithID, or any other type.
     * @param id - numerical id (Integer).
     */
    constructor(valueOrObject: ValueOrObject<T>, id?: number);
    /**
     * Serialize class instance to JSON.
     * @example {"id" : 0, "value" : "Si" }
     */
    toJSON(): ObjectWithIdAndValue<T>;
}
