import { ArrayWithIds } from "@mat3ra/code";
import { AtomicConstraintSchema } from "@mat3ra/esse/dist/js/types";
export type AtomicConstraintValue = AtomicConstraintSchema["value"];
export declare class AtomicConstraints extends ArrayWithIds<AtomicConstraintValue> {
    /**
     * Get constraints for an atom with index as string.
     * @param idx - atom index.
     * @param mapFn (OPTIONAL) - a function to be applied to each constraint. By default 0 or 1 is returned.
     */
    getAsStringByIndex(idx: number, mapFn?: (val: boolean) => string): string;
}
