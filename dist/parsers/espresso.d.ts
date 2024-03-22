import { MaterialSchema } from "@mat3ra/esse/lib/js/types";
/**
 * Construct textual representation of a materialOrConfig according to Quantum ESPRESSO pw.x input format.
 * @param materialOrConfig - material class instance or its config object
 */
declare function toEspressoFormat(materialOrConfig: MaterialSchema): string;
declare const _default: {
    toEspressoFormat: typeof toEspressoFormat;
};
export default _default;
