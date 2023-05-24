import s from "underscore.string";

import { ConstrainedBasis } from "../basis/constrained_basis";
import { ATOMIC_COORD_UNITS, units } from "../constants";
import { Lattice } from "../lattice/lattice";
import math from "../math";
import { Material } from "../material";
import { LATTICE_TYPE } from "../lattice/types";

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
export function atomsCount(poscarFileContent) {
    const atomsLine = poscarFileContent.split("\n")[6].split(" ");
    return atomsLine.map((x) => parseInt(x, 10)).reduce((a, b) => a + b);
}

/**
 * Parses POSCAR file into a Material object.
 * @param {string} fileContent - POSCAR file content.
 * @return {Object} Material Config.
 */
function fromPoscar(fileContent) {
    const lines = fileContent.split("\n");

    const comment = lines[0];
    const latticeConstant = parseFloat(lines[1].trim());

    const latticeVectors = [
        lines[2].split(" ").map(Number),
        lines[3].split(" ").map(Number),
        lines[4].split(" ").map(Number),
    ].map((vector) => vector.map((component) => component * latticeConstant));

    // Atom symbols and counts
    const atomSymbols = lines[5].trim().split(" ");
    const atomCounts = lines[6].trim().split(" ").map(Number);

    // Check if selective dynamics and coordinates type is used
    let selectiveDynamics = false;
    let coordinateType = "";
    let startLine = 7;
    if (lines[startLine].trim()[0].toLowerCase() === "s") {
        selectiveDynamics = true;
        coordinateType = lines[8].trim().toLowerCase();
        startLine = 9;
    } else {
        coordinateType = lines[7].trim().toLowerCase();
        startLine = 8;
    }

    // Atom coordinates and constraints
    const coordinates = [];
    const constraints = [];
    let atomIndex = 0;
    for (let i = 0; i < atomSymbols.length; i++) {
        for (let j = 0; j < atomCounts[i]; j++) {
            const lineComponents = lines[startLine + atomIndex].trim().split(" ");
            const coordinate = [
                parseFloat(lineComponents[0]),
                parseFloat(lineComponents[1]),
                parseFloat(lineComponents[2]),
            ];
            coordinates.push(coordinate);
            if (selectiveDynamics) {
                const constraint = [
                    lineComponents[3] === "T",
                    lineComponents[4] === "T",
                    lineComponents[5] === "T",
                ];
                constraints.push(constraint);
            }
            atomIndex++;
        }
    }
    
    const latticeType = LATTICE_TYPE.HEX; // TODO: handle actual lattice

    const materialConfig = {
        name: comment,
        basis: {
            elements: atomSymbols.flatMap((symbol, i) =>
                Array(atomCounts[i]).fill(symbol)
            ),
            coordinates,
            units: coordinateType === "direct" ? ATOMIC_COORD_UNITS.direct : ATOMIC_COORD_UNITS.cartesian,
            cell:  latticeVectors,
        },
        lattice: new Lattice(latticeType, latticeVectors)
    };

    if (selectiveDynamics) {
        materialConfig.basis.constraints = constraints;
    }

    return materialConfig;
}


export default {
    toPoscar,
    fromPoscar,
    atomicConstraintsCharFromBool,
    atomsCount,
};
