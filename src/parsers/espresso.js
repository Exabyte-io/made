import _ from "underscore";
import s from "underscore.string";

import { ConstrainedBasis } from "../basis/constrained_basis";
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
            result.card_lines = lines;
        }
    });
    return result;
}

/**
 * @summary Function to get atomic species from Fortran namelist.
 * @param {String[]} lines - Fortran card_lines string.
 * @param {Number} nSpecies - Number of atomic species.
 * @return {String[]} - atomic species labels.
 */
function getAtomicSpecies(lines, nSpecies) {
    // filter out blank or comment lines
    const trimmedLines = lines.filter((line) => line.trim() && !line.startsWith("#"));
    const species = [];

    trimmedLines.forEach((line) => {
        // TODO: redo here everything. Something is uttterly wrong
        if (line.startsWith("ATOMIC_SPECIES")) {
            // find the next ATOMIC_SPECIES line
            const atomicSpeciesLineIndex = trimmedLines.findIndex((_line) =>
                _line.startsWith("ATOMIC_SPECIES"),
            );

            for (let i = 0; i < nSpecies; i++) {
                const labelWeightPseudo = trimmedLines[atomicSpeciesLineIndex + i + 1].split(" ");
                species.push({
                    label: labelWeightPseudo[0],
                    // Only atom symbol is used now
                    // weight: parseFloat(labelWeightPseudo[1]),
                    // pseudo: labelWeightPseudo[2],
                });
            }
        }
    });
    return species;
}

/**
 * @summary Parse unit cell from CELL_PARAMETERS card.
 * @param {String[]} lines - Fortran card_lines string.
 * @param {Number | Null} alat -
 * @returns {Number[][]} cell - unit cell vectors.
 */
function getCellParameters(lines, alat = null) {
    let cell = null;
    let cellAlat = null;

    const trimmedLines = lines.filter((line) => line.trim() && line[0] !== "#");
    let cellUnits;

    trimmedLines.forEach((line, index) => {
        if (line.trim().startsWith("CELL_PARAMETERS")) {
            if (cell !== null) {
                throw new Error("CELL_PARAMETERS specified multiple times");
            }
            if (line.toLowerCase().includes("bohr")) {
                if (alat !== null) {
                    throw new Error(
                        "Lattice parameters given in &SYSTEM celldm/A and CELL_PARAMETERS bohr",
                    );
                }
                cellUnits = coefficients.BOHR_TO_ANGSTROM;
            } else if (line.toLowerCase().includes("angstrom")) {
                if (alat !== null) {
                    throw new Error(
                        "Lattice parameters given in &SYSTEM celldm/A and CELL_PARAMETERS angstrom",
                    );
                }
                cellUnits = 1.0;
            } else if (line.toLowerCase().includes("alat")) {
                if (line.includes("=")) {
                    // eslint-disable-next-line no-param-reassign
                    alat = parseFloat(line.trim().split("=")[1]) * coefficients.BOHR_TO_ANGSTROM;
                    cellAlat = alat;
                } else if (alat === null) {
                    throw new Error("Lattice parameters must be set in &SYSTEM for alat units");
                }
                cellUnits = alat;
            } else if (alat === null) {
                cellUnits = coefficients.BOHR_TO_ANGSTROM;
            } else {
                cellUnits = alat;
            }
            // Grab the parameters; blank lines have been removed
            cell = [];
            for (let i = 0; i < 3; i++) {
                cell.push(trimmedLines[index + 1 + i].split(/\s+/).slice(0, 3).map(parseFloat));
            }
            cell = cell.map((row) => row.map((value) => value * cellUnits));
        }
    });
    return [cell, cellAlat];
}

/**
 * Parses QE .in file into a Material config object.
 * @param {String} fileContent - contents of the .in file
 * @return {Object} Material config.
 */
function fromEspressoFormat(fileContent) {
    const data = parseFortranNameList(fileContent);
    const nSpecies = data.system.nat;
    const cell = getCellParameters(data.card_lines);
    const elements = getAtomicSpecies(data.card_lines, nSpecies);

    const latticeVectors = cell;
    const coordinates = [];

    const lattice = Lattice.fromVectors(latticeVectors);

    const basis = new ConstrainedBasis({
        elements,
        coordinates,
        units: ATOMIC_COORD_UNITS.crystal,
        cell: lattice.vectorArrays,
    });

    const materialConfig = {
        lattice: lattice.toJSON(),
        basis: basis.toJSON(),
        name: "comment", // TODO: add actual system name
        isNonPeriodic: false,
    };
    return materialConfig;
}

export default {
    isEspressoFormat,
    toEspressoFormat,
    fromEspressoFormat,
};
