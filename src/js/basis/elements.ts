import { ArrayWithIds } from "@mat3ra/code";
import { AtomicElementSchema } from "@mat3ra/esse/dist/js/types";
import { uniq } from "lodash";

export type AtomicElementValue = AtomicElementSchema["value"];

export class Elements extends ArrayWithIds<AtomicElementValue> {
    getUnique(): string[] {
        return uniq(this.values.map((element) => element)) as string[];
    }
}
