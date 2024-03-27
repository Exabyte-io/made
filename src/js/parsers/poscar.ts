import s from "underscore.string";

import { ConstrainedBasis } from "../basis/constrained_basis";
import { Coordinate } from "../basis/types";
import { ATOMIC_COORD_UNITS } from "../constants";
import { Constraint, ConstraintValue } from "../constraints/constraints";
import { Lattice } from "../lattice/lattice";
import { Vector } from "../lattice/types";
import math from "../math";
import { MaterialJSON } from "../types";

const _print = (x: number, printFormat = "%14.9f") => s.sprintf(printFormat, math.precise(x));
const _latticeVectorsToString = (vectors: Vector[]) =>
    vectors.map((v) => v.map((c) => _print(c)).join("\t")).join("\n");
const atomicConstraintsCharFromBool = (bool: boolean): string => (bool ? "T" : "F");

/**
 * Obtain a textual representation of a material in POSCAR format.
 * @param materialOrConfig - material class instance or config object.
 * @param omitConstraints - whether to discard constraints passed with material.
 */
function toPoscar(materialOrConfig: MaterialJSON, omitConstraints = false): string {
    const lattice = new Lattice(materialOrConfig.lattice);
    const vectorsAsString = _latticeVectorsToString(lattice.vectorArrays);
    // @ts-ignore
    const basis = new ConstrainedBasis({
        ...materialOrConfig.basis,
        cell: lattice.vectorArrays,
    });
    const BasisLines: string[] = [];
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
    const countsLine = basis.elementCounts.map((e) => parseInt(`${e.count}`, 10)).join(" ");
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
 */
export function atomsCount(poscarFileContent: string): number {
    const atomsLine = poscarFileContent.split("\n")[6].split(/\s+/);
    return atomsLine.map((x) => parseInt(x, 10)).reduce((a, b) => a + b);
}

/**
 * Parses POSCAR file into a Material config object.
 * @param fileContent - POSCAR file content.
 * @return Material config.
 */
function fromPoscar(fileContent: string): object {
    const lines = fileContent.split("\n");

    const comment = lines[0];
    const latticeConstant = parseFloat(lines[1].trim());

    // Atom symbols and counts
    const atomSymbols = lines[5].trim().split(/\s+/);
    const atomCounts = lines[6].trim().split(/\s+/).map(Number);

    // Check if selective dynamics and coordinates type is used
    let selectiveDynamics = false;
    let coordinateType = "";
    let startLine = 7;
    if (lines[startLine].trim()[0].toLowerCase() === "s") {
        selectiveDynamics = true;
        coordinateType = lines[8].trim()[0].toLowerCase();
        startLine = 9;
    } else {
        coordinateType = lines[7].trim()[0].toLowerCase();
        startLine = 8;
    }

    const elements: string[] = atomSymbols
        .map((symbol: string, i: number): string[] => Array(atomCounts[i]).fill(symbol))
        .reduce((a, b) => a.concat(b), []);

    // Atom coordinates and constraints
    const coordinates: Coordinate[] = [];
    const constraints: Constraint[] = [];
    let atomIndex = 0;
    for (let i = 0; i < atomSymbols.length; i++) {
        for (let j = 0; j < atomCounts[i]; j++) {
            const lineComponents = lines[startLine + atomIndex].trim().split(/\s+/);
            const coordinate: Coordinate = [
                parseFloat(lineComponents[0]),
                parseFloat(lineComponents[1]),
                parseFloat(lineComponents[2]),
            ];
            coordinates.push(coordinate);

            // Add constraints if selective dynamics is used
            if (selectiveDynamics) {
                const constraint: ConstraintValue = [
                    lineComponents[3] === "T",
                    lineComponents[4] === "T",
                    lineComponents[5] === "T",
                ];
                constraints.push({ id: j, value: constraint });
            }
            atomIndex += 1;
        }
    }

    const lattice = Lattice.fromVectors({
        a: lines[2].trim().split(/\s+/).map(Number) as Vector,
        b: lines[3].trim().split(/\s+/).map(Number) as Vector,
        c: lines[4].trim().split(/\s+/).map(Number) as Vector,
        alat: latticeConstant,
        units: "angstrom",
    });

    const basis = new ConstrainedBasis({
        elements,
        coordinates,
        units: coordinateType === "c" ? ATOMIC_COORD_UNITS.cartesian : ATOMIC_COORD_UNITS.crystal,
        cell: lattice.vectorArrays,
        constraints,
    });

    const materialConfig = {
        lattice: lattice.toJSON(),
        basis: basis.toJSON(),
        name: comment,
        isNonPeriodic: false,
    };

    return materialConfig;
}

/**
 * @summary Checks if a string has a POSCAR format (first 8 lines are read)
 * @param text - string to check
 */
function isPoscar(text: string): boolean {
    const lines = text.split("\n");

    // Checking number of lines, minimum requirement for POSCAR
    if (lines.length < 7) {
        return false;
    }

    // Check the lattice constant (a single floating point number)
    if (!/^[+-]?[\d.]+$/.test(lines[1].trim())) {
        return false;
    }

    // Check the lattice vectors (three lines, each with three floating point numbers)
    for (let i = 2; i <= 4; i++) {
        if (!/^[-+]?\d*\.\d+\s+[-+]?\d*\.\d+\s+[-+]?\d*\.\d+$/.test(lines[i].trim())) {
            return false;
        }
    }

    // Check the atomic species line (alphabetic characters, space-separated)
    if (!/^[a-zA-Z\s]+$/.test(lines[5].trim())) {
        return false;
    }

    // Check the number of atoms per species line (digits, space-separated)
    if (!/^[\d\s]+$/.test(lines[6].trim())) {
        return false;
    }

    // Check the coordinate type line (only first character is neccessary, "s" or "S" for "Selective dynamics",
    // "d" or "D" for "Direct", or "c" or "C" for "Cartesian")
    if (!/^[sdc]/.test(lines[7].trim().toLowerCase())) {
        return false;
    }
    return true;
}

export default {
    isPoscar,
    toPoscar,
    fromPoscar,
    atomicConstraintsCharFromBool,
    atomsCount,
};
