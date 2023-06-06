import _ from "underscore";
import s from "underscore.string";

import { ConstrainedBasis } from "../basis/constrained_basis";
// eslint-disable-next-line no-unused-vars
import { ATOMIC_COORD_UNITS, coefficients } from "../constants";
import { Lattice } from "../lattice/lattice";
import xyz from "./xyz";

/**
 * Construct textual representation of a materialOrConfig according to Quantum ESPRESSO pw.x input format.
 * @param {Material|Object} materialOrConfig - material class instance or its config object
 * @return {String}
 */
function toEspressoFormat(materialOrConfig) {
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

/**
 * @summary checks if the given fileContent is in the format of a Quantum ESPRESSO .in file.
 * @param {String} fileContent
 * @return {boolean}
 */
function isEspressoFormat(fileContent) {
    const lines = fileContent.split("\n");
    if (lines[0].trim() === "&CONTROL") return true;
    return false;
}

/**
 * @summary Function to convert Fortran namelist to JavaScript object.
 * @param {String} data - Fortran string.
 * @return {Object}
 */
function parseFortranNameList(data) {
    const blocks = data.split(/\/\n+/);
    const result = {};

    blocks.forEach((block) => {
        const lines = block.split("\n");
        const blockName = lines[0].trim().toLowerCase();

        if (blockName.startsWith("&")) {
            const blockContent = {};
            for (let i = 1; i < lines.length; i++) {
                const line = lines[i].trim();
                if (line) {
                    const [key, value] = line.split("=").map((str) => str.trim());
                    // convert Fortran boolean to JavaScript boolean
                    if (value === ".true.") {
                        blockContent[key] = true;
                    } else if (value === ".false.") {
                        blockContent[key] = false;
                    } else if (Number.isNaN(Number(value))) {
                        // if value is not a number, keep it as a string
                        blockContent[key] = value.replace(/'/g, ""); // remove single quotes
                    } else {
                        blockContent[key] = Number(value); // convert to number
                    }
                }
            }
            result[blockName.slice(1)] = blockContent; // remove '&' from block name
        } else {
            result.cards = block;
        }
    });
    return result;
}

// function parseCards(cards) {
//     const result = {};
//     const sections = cards.split(
//         /^\s*(ATOMIC_SPECIES|ATOMIC_POSITIONS|K_POINTS|CELL_PARAMETERS|OCCUPATIONS|CONSTRAINTS)\s/gm,
//     );
//
//     sections.forEach((section) => {
//         const sectionName = section.trim();
//         const sectionData = section
//             .trim()
//             .split("\n")
//             .map((line) => {
//                 // Remove any comments and trim line
//                 const cleanLine = line.replace(/(!|#).*$/, "").trim();
//                 if (cleanLine) {
//                     // split by spaces
//                     return cleanLine.split(/\s+/);
//                 }
//                 return null;
//             })
//             .filter(Boolean); // Remove null entries
//
//         result[sectionName.convertKeysToCamelCaseForObject()] = sectionData;
//     });
//
//     return result;
// }

// /**
//  * @summary Function to get atomic species from Fortran namelist.
//  * @param {String} cards - Fortran cards string.
//  * @param {Number} nSpecies - Number of atomic species.
//  * @return {{symbol: String, weight: Number, pseudo: String}[]} - atomic species labels.
//  */
// function getAtomicSpecies(cards, nSpecies) {
//     const species = [];
//     const regex = /ATOMIC_SPECIES\s+((?:\S+\s+\S+\s+\S+\s*){1,nSpecies})/s;
//     const match = cards.match(regex);
//     console.log(match);
//
//     if (match && match[1]) {
//         const speciesLines = match[1].trim().split("\n");
//
//         for (const line of speciesLines) {
//             const [label, weight, pseudo] = line.trim().split(/\s+/);
//             species.push({
//                 label,
//                 weight: parseFloat(weight),
//                 pseudo,
//             });
//         }
//     }
//
//     return species;
// }
//
// /**
//  * @summary Parse unit cell from CELL_PARAMETERS card.
//  * @param {String} cards - Fortran card_lines string.
//  * @param {Number | Null} alat
//  * @return {{Number[][], Number}} [cell, cellAlat] - unit cell vectors.
//  */
// function getCellParameters(cards, alat = null) {
//     // Regular expression to find cell parameters section
//     // const regex = /CELL_PARAMETERS\s*\(([^)]+)\)((?:[\s\S](?!#))*)/s;
//     // const match = cards.match(regex);
//     // console.log(match);
//     // if (!match) {
//     //     throw new Error("CELL_PARAMETERS not found in input");
//     // }
//
//     const regex =
//         /ATOMIC_POSITIONS \(w*?\)([\s\S]*?)CELL_PARAMETERS \(w*?\)([\s\S]*?)K_POINTS automatic([\s\S]*)/;
//
//     const match = input.match(regex);
//
//     if (match) {
//         const [, positions, parameters, kPoints] = match;
//
//         const atomicSpecies = {
//             positions: positions
//                 .trim()
//                 .split("\n")
//                 .map((line) => line.trim().split(" ").slice(1).map(Number)),
//             parameters: parameters
//                 .trim()
//                 .split("\n")
//                 .map((line) => line.trim().split(" ").map(Number)),
//             kPoints: kPoints.trim().split(" ").map(Number),
//         };
//
//         console.log(atomicSpecies);
//     } else {
//         console.log("No match found.");
//     }
//
//     let cellUnits, cellAlat;
//     const units = match[1].toLowerCase();
//
//     if (units.includes("bohr")) {
//         if (alat !== null) {
//             throw new Error(
//                 "Lattice parameters given in &SYSTEM celldm/A and CELL_PARAMETERS bohr",
//             );
//         }
//         cellUnits = coefficients.BOHR_TO_ANGSTROM;
//     } else if (units.includes("angstrom")) {
//         if (alat !== null) {
//             throw new Error(
//                 "Lattice parameters given in &SYSTEM celldm/A and CELL_PARAMETERS angstrom",
//             );
//         }
//         cellUnits = 1.0;
//     } else if (units.includes("alat")) {
//         if (match[1].includes("=")) {
//             cellAlat = parseFloat(units.split("=")[1]) * coefficients.BOHR_TO_ANGSTROM;
//         } else if (alat === null) {
//             throw new Error("Lattice parameters must be set in &SYSTEM for alat units");
//         }
//         cellUnits = alat;
//     } else if (alat === null) {
//         cellUnits = coefficients.BOHR_TO_ANGSTROM;
//     } else {
//         cellUnits = alat;
//     }
//
//     const lines = match[2]
//         .split("\n")
//         .filter((line) => line.trim() && !line.trim().startsWith("#"));
//
//     if (lines.length < 3) {
//         throw new Error("Incomplete CELL_PARAMETERS data");
//     }
//
//     const cell = lines.slice(0, 3).map((str) => {
//         return str
//             .split(/\s+/)
//             .filter((el) => el)
//             .map((num) => parseFloat(num) * cellUnits);
//     });
//
//     return { cell, cellAlat };
// }
//
// /**
//  * @summary Parses atomic positions from ATOMIC_POSITIONS card.
//  * @param {String[]} lines
//  * @param {Number} nAtoms
//  * @param {Number[][] | null} cell - unit cell with angstroms
//  * @param {Number | null} alat
//  * @return {{coordinates: Number[][], constraints: Object[]}}
//  */
// function getAtomicPositions(lines, nAtoms, cell = null, alat = null) {
//     const coordinates = [];
//     const constraints = [];
//
//     const trimmedLines = lines.filter((line) => line.trim() && line[0] !== "#");
//
//     trimmedLines.forEach((line, index) => {
//         if (line.trim().startsWith("ATOMIC_POSITIONS")) {
//             for (let i = 0; i < nAtoms; i++) {
//                 const [element, x, y, z, ifX, ifY, ifZ] = trimmedLines[index + 1 + i].split(/\s+/);
//
//                 coordinates.push([parseFloat(x), parseFloat(y), parseFloat(z)]);
//                 constraints.push({ id: i, value: [ifX === "T", ifY === "T", ifZ === "T"] });
//             }
//         }
//     });
//
//     return { coordinates, constraints };
// }

/**
 * @summary Read unit cell parameters from CELL_PARAMETERS card
 * @param {String} cards
 * @param {String} cell
 * @returns {{Number[][], Number}}
 */

// eslint-disable-next-line no-unused-vars
function getCellParameters(cards, cell) {
    const cellParamsRegex = /CELL_PARAMETERS\s\(([\w]+)\)(?:\n?([!#].*)?)\n(((-?\d+.\d+)\s+)*)/gm;
    // Extract cell parameters section
    const cellParamsMatch = cellParamsRegex.exec(cards);

    console.log(cellParamsMatch);
    // Extract the unit
    // eslint-disable-next-line no-unused-vars
    const unit = cellParamsMatch[1];

    const cellParamsData = cellParamsMatch[3].trim().split("\n");

    const cellVectors = [];
    cellParamsData.forEach((line) => {
        const values = line
            .trim()
            .match(/([\d.-]+)/g)
            .map(Number);
        cellVectors.push(values);
    });
    console.log("cellvectors:::::", cellVectors);

    const alat = 1.0; // TODO: depends on unit
    return {
        cell: cellVectors,
        alat,
    };
}

/**
 * @summary function to get atomic positions from the ATOMIC_POSITIONS card
 * @param {String} cards
 * @param {Number} nAtoms
 * @returns {[String], [Number], [String], [Boolean]}
 */
function getAtomicPositions(cards) {
    const regex =
        /ATOMIC_POSITIONS\s\(([\w]+)\)(?:\n?([!#].*)?)\n((?:\w+\s+[-?\d+.\d+\s+]*\n)*(?!\w+\s+\())?/gm;
    const match = regex.exec(cards);

    if (!match) {
        throw new Error("ATOMIC_POSITIONS section not found in input data.");
    }

    // const coordSystem = match[1]; // TODO: use that info to multiply appropriayle
    const atomicPosData = match[3];

    const elements = [];
    const coordinates = [];
    const constraints = [];

    atomicPosData.split("\n").forEach((line) => {
        if (line.trim()) {
            const [atom, x, y, z, if_x, if_y, if_z] = line.trim().split(/\s+/);
            elements.push(atom);
            coordinates.push([Number(x), Number(y), Number(z)]);
            constraints.push([Boolean(Number(if_x)), Boolean(Number(if_y)), Boolean(Number(if_z))]);
        }
    });

    return {
        elements,
        coordinates,
        constraints,
    };
}

/**
 * Parses QE .in file into a Material config object.
 * @param {String} fileContent - contents of the .in file
 * @return {Object} Material config.
 */
function fromEspressoFormat(fileContent) {
    const data = parseFortranNameList(fileContent);
    // const nSpecies = data.system.ntyp;
    // const nAtoms = data.system.nat;
    // const ibrav = data.system.ibrav;
    const { cards } = data;
    console.log(cards);

    const cell = getCellParameters(cards);
    const { elements } = getAtomicPositions(cards);
    const { coordinates } = getAtomicPositions(cards);
    const { constraints } = getAtomicPositions(cards);

    const lattice = Lattice.fromVectors({
        a: cell.cell[0],
        b: cell.cell[1],
        c: cell.cell[2],
        alat: cell.alat,
        units: ATOMIC_COORD_UNITS.angstrom,
    });

    const basis = new ConstrainedBasis({
        elements,
        coordinates,
        units: ATOMIC_COORD_UNITS.cartesian,
        cell: lattice.vectorArrays,
        constraints,
    });

    const materialConfig = {
        lattice: lattice.toJSON(),
        basis: basis.toJSON(),
        name: data.control.title,
        isNonPeriodic: false,
    };
    return materialConfig;
}

export default {
    isEspressoFormat,
    toEspressoFormat,
    fromEspressoFormat,
};
