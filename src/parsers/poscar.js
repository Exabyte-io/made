import s from "underscore.string";

import { ConstrainedBasis } from "../basis/constrained_basis";
import { ATOMIC_COORD_UNITS } from "../constants";
import { Lattice } from "../lattice/lattice";
import math from "../math";

const _print = (x, printFormat = "%14.9f") => s.sprintf(printFormat, math.precise(x));

const _latticeVectorsToString = (vectors) =>
    vectors.map((v) => v.map((c) => _print(c)).join("\t")).join("\n");

const atomicConstraintsCharFromBool = (bool) => (bool ? "T" : "F");

/**
 * Obtain a textual representation of a material in POSCAR format.
 * @param {Material|Object} materialOrConfig - material class instance or config object.
 * @param {Boolean} omitConstraints - whether to discard constraints passed with material.
 * @return {string}
 */
function toPoscar(materialOrConfig, omitConstraints = false) {
    const lattice = new Lattice(materialOrConfig.lattice);
    const vectorsAsString = _latticeVectorsToString(lattice.vectorArrays);

    const basis = new ConstrainedBasis({
        ...materialOrConfig.basis,
        cell: lattice.vectorArrays,
    });

    const BasisLines = [];
    let addSelectiveDynamics = false;
    basis._elements.array.forEach((item, idx) => {
        const coord = basis.getCoordinateByIndex(idx).map((x) => _print(x));
        const constraintsAsString = omitConstraints
            ? ""
            : basis.AtomicConstraints.getAsStringByIndex(idx, atomicConstraintsCharFromBool);
        if (constraintsAsString && !omitConstraints) addSelectiveDynamics = true;
        BasisLines.push([coord.join(" "), constraintsAsString, item].join(" "));
    });
    const basisContent = BasisLines.join("\n");

    const elementsLine = basis.elementCounts.map((e) => e.value).join(" ");
    const countsLine = basis.elementCounts.map((e) => parseInt(e.count, 10)).join(" ");

    const coordsType =
        materialOrConfig.basis.units === ATOMIC_COORD_UNITS.cartesian ? "cartesian" : "direct";

    return [
        materialOrConfig.name,
        "1.0",
        vectorsAsString,
        elementsLine,
        countsLine,
        // add selective dynamics only if there are some constraints!
        ...(addSelectiveDynamics ? ["Selective dynamics"] : []),
        coordsType,
        basisContent,
    ].join("\n");
}

export default {
    toPoscar,
    atomicConstraintsCharFromBool,
};
