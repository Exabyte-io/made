import { ArrayWithIds } from "@mat3ra/code";
import { AtomicElementSchema } from "@mat3ra/esse/dist/js/types";
import * as _ from "underscore";

export type AtomicElementValue = AtomicElementSchema["value"];

export class Elements extends ArrayWithIds<AtomicElementValue> {
    getUnique(): string[] {
        return _.unique(this.values.map((element) => element)) as string[];
    }
}
