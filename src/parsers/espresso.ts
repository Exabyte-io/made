import { MaterialSchema } from "@mat3ra/esse/lib/js/types";
import _ from "underscore";
import s from "underscore.string";

import { Lattice } from "../lattice/lattice";
import xyz from "./xyz";

/**
 * Construct textual representation of a materialOrConfig according to Quantum ESPRESSO pw.x input format.
 * @param materialOrConfig - material class instance or its config object
 */
function toEspressoFormat(materialOrConfig: MaterialSchema): string {
    const l = new Lattice(materialOrConfig.lattice);
    const vectors = l.vectorArrays;
    const vectorsAsString = _.map(vectors, (v) => {
        return `${s.sprintf("%14.9f", v[0])}\t${s.sprintf("%14.9f", v[1])}\t${s.sprintf(
            "%14.9f",
            v[2],
        )}`;
    }).join("\n");
    return s.sprintf(
        "CELL_PARAMETERS (angstroms)\n%s\n\nATOMIC_POSITIONS (crystal)\n%s",
        vectorsAsString,
        xyz.fromMaterial(materialOrConfig),
    );
}

export default {
    toEspressoFormat,
};
