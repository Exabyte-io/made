import { ArrayWithIds, ValueWithId } from "@mat3ra/code";
import { AtomicElementSchema } from "@mat3ra/esse/dist/js/types";
export type AtomicElementValue = AtomicElementSchema["value"];
export declare class Element extends ValueWithId<AtomicElementValue> {
    value: AtomicElementValue;
    constructor({ value, id }: AtomicElementSchema);
    prettyPrint(): string;
}
export declare class Elements extends ArrayWithIds<AtomicElementValue> {
    getUnique(): string[];
}
