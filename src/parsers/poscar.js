import s from "underscore.string";

import { ConstrainedBasis } from "../basis/constrained_basis";
import { ATOMIC_COORD_UNITS } from "../constants";
import { defaultMaterialConfig } from "../default_material";
import { Lattice } from "../lattice/lattice";
import math from "../math";
import xyz from "./xyz";

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

/**
 * @summary calculates the number of atoms in a poscar file based on the summation of the numbers in line 7 of the file.
 * Poscar file formatting: https://www.vasp.at/wiki/index.php/POSCAR
 * @param poscarFileContent
 * @returns {Number}
 */
export function getAtomsCount(poscarFileContent) {
    const atomsLine = poscarFileContent.split("\n")[6].split(" ");
    return atomsLine.map((x) => parseInt(x, 10)).reduce((a, b) => a + b);
}

/**
 * Converts XYZ formated string content to POSCAR formated string content.
 * @param xyzContent
 * @returns {string}
 */
export function convertFromXyz(xyzContent) {
    const xyzConfig = defaultMaterialConfig;
    const xyzArray = xyzContent.split(/\r?\n/);
    const xyzArrayBasisOnly = xyzArray.slice(2, -1);
    const xyzBasis = xyzArrayBasisOnly.join("\n");
    xyzConfig.basis = xyz.toBasisConfig(xyzBasis);
    xyzConfig.basis.units = "cartesian";
    return toPoscar(xyzConfig);
}

/**
 * Function converts a string of content into a POSCAR formatted string of content. The function called to convert the
 * format is determined by the format type passed to the function.
 * @param {String} xyzContent
 * @param {String} format
 * @returns {String}
 */
export function convertFromOtherFormat(content, format) {
    if (format === "xyz") {
        return convertFromXyz(content);
    }
    return content;
}

export default {
    toPoscar,
    atomicConstraintsCharFromBool,
    getAtomsCount,
    convertFromOtherFormat,
    convertFromXyz,
};
