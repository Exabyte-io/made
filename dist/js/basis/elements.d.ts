import { ArrayWithIds } from "@mat3ra/code";
import { AtomicElementSchema } from "@mat3ra/esse/dist/js/types";
export type AtomicElementValue = AtomicElementSchema["value"];
export declare class Elements extends ArrayWithIds<AtomicElementValue> {
    getUnique(): string[];
}
