import { ArrayWithIds } from "@mat3ra/code";
import { ElementSchema } from "@mat3ra/esse/dist/js/types";
import * as _ from "underscore";

/**
 * Specialized class for handling atomic elements in a basis
 * Extends ArrayWithIds to provide element-specific functionality
 */
export class Elements extends ArrayWithIds<string> {
    getUniqueElements(): string[] {
        return _.unique(this.values);
    }
}
