import { ArrayWithIds } from "@mat3ra/code";
import * as _ from "underscore";

export class Elements extends ArrayWithIds<string> {
    getUnique(): string[] {
        return _.unique(this.values.map((element) => element)) as string[];
    }
}
