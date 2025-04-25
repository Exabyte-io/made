import { ArrayWithIds, ValueWithId } from "@mat3ra/code";
import { AtomicConstraintSchema } from "@mat3ra/esse/dist/js/types";
import { padEnd } from "lodash";

export type AtomicConstraintValue = AtomicConstraintSchema["value"];

export class Constraint extends ValueWithId<AtomicConstraintValue> {
    value: AtomicConstraintValue;

    constructor({ value, id }: AtomicConstraintSchema) {
        super({ id, value });
        this.value = value;
    }

    getValueAsString(): string {
        return this.value.map((x) => (x ? 1 : 0)).join(" ");
    }

    prettyPrint(
        element: string,
        coordinate: { prettyPrint: () => string },
        constraint: boolean[],
    ): string {
        return (
            padEnd(element, 4) +
            coordinate.prettyPrint() +
            " " +
            constraint.map((x) => (x ? 1 : 0)).join(" ")
        );
    }
}

export class AtomicConstraints extends ArrayWithIds<AtomicConstraintValue> {
    /**
     * Get constraints for an atom with index as string.
     * @param idx - atom index.
     * @param mapFn (OPTIONAL) - a function to be applied to each constraint. By default 0 or 1 is returned.
     */
    getAsStringByIndex(idx: number, mapFn = (val: boolean): string => (val ? "1" : "0")): string {
        const constraints = this.getElementValueByIndex(idx);
        return constraints ? constraints.map(mapFn).join(" ") : "";
    }
}
