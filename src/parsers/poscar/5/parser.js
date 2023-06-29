import { ATOMIC_COORD_UNITS } from "@exabyte-io/code.js/dist/constants";

import { ConstrainedBasis } from "../../../basis/constrained_basis";
import { Lattice } from "../../../lattice/lattice";

/**
 * @summary Checks if a string has a POSCAR format (first 8 lines are read)
 * @param {string} text - string to check
 * @returns {boolean}
 */
function isPoscar(text) {
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
        if (!/^[+-]?[\d.\s]+$/.test(lines[i].trim())) {
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

/**
 * Parses POSCAR file into a Material config object.
 * @param {string} fileContent - POSCAR file content.
 * @return {Object} Material config.
 */
function fromPoscar(fileContent) {
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

    const elements = atomSymbols
        .map((symbol, i) => Array(atomCounts[i]).fill(symbol))
        .reduce((a, b) => a.concat(b), []);

    // Atom coordinates and constraints
    const coordinates = [];
    const constraints = [];
    let atomIndex = 0;
    for (let i = 0; i < atomSymbols.length; i++) {
        for (let j = 0; j < atomCounts[i]; j++) {
            const lineComponents = lines[startLine + atomIndex].trim().split(/\s+/);
            const coordinate = [
                parseFloat(lineComponents[0]),
                parseFloat(lineComponents[1]),
                parseFloat(lineComponents[2]),
            ];
            coordinates.push(coordinate);

            // Add constraints if selective dynamics is used
            if (selectiveDynamics) {
                const constraint = [
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
        a: lines[2].trim().split(/\s+/).map(Number),
        b: lines[3].trim().split(/\s+/).map(Number),
        c: lines[4].trim().split(/\s+/).map(Number),
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

    const cell = lattice;
    const { units } = basis;
    const name = comment;
    console.log(cell, elements, coordinates, constraints, units, name);
    return { cell, elements, coordinates, constraints, units, name };
}

export default {
    fromPoscar,
    isPoscar,
};
