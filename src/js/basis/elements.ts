import { ArrayWithIds, ValueWithId } from "@mat3ra/code";
import { AtomicElementSchema } from "@mat3ra/esse/dist/js/types";
import { uniq } from "lodash";

export type AtomicElementValue = AtomicElementSchema["value"];

export class Element extends ValueWithId<AtomicElementValue> {
    value: AtomicElementValue;

    constructor({ value, id }: AtomicElementSchema) {
        super({ id, value });
        this.value = value;
    }

    prettyPrint(): string {
        return this.value;
    }
}

export class Elements extends ArrayWithIds<AtomicElementValue> {
    getUnique(): string[] {
        return uniq(this.values.map((element) => element)) as string[];
    }
}
