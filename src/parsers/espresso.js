import _ from "underscore";
import s from "underscore.string";

import { ConstrainedBasis } from "../basis/constrained_basis";
// eslint-disable-next-line no-unused-vars
import { ATOMIC_COORD_UNITS, coefficients } from "../constants";
import { Lattice } from "../lattice/lattice";
import math from "../math";
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
                    // Check if the line contains celldm(i)
                    const celldmMatch = line.match(/celldm\((\d+)\)\s+=\s+(.+)/);
                    if (celldmMatch) {
                        const index = Number(celldmMatch[1]);
                        const value = Number(celldmMatch[2]);
                        if (!blockContent.celldm) {
                            blockContent.celldm = {}; // Initialize if not present
                        }
                        blockContent.celldm[index] = value; // Store value at index
                    } else {
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
            }
            result[blockName.slice(1)] = blockContent; // remove '&' from block name
        } else {
            result.cards = block;
        }
    });

    return result;
}

/**
 * @summary Read unit cell parameters from CELL_PARAMETERS card
 * @param {String} cards
 * @param {Number} alat
 * @returns {cell: Number[][], alat: Number, units: String}
 */
function getCellParameters(cards, alat = null) {
    const cellParamsRegex = /CELL_PARAMETERS\s\(([\w]+)\)(?:\n?([!#].*)?)\n(((-?\d+.\d+)\s+)*)/gm;
    // Extract cell parameters section
    const cellParamsMatch = cellParamsRegex.exec(cards);

    // Extract the units
    const units = cellParamsMatch[1];
    let cellUnits = 1.0;

    if (units === "bohr") {
        if (alat !== null) {
            throw new Error(
                "Lattice parameters given in &SYSTEM celldm/A and CELL_PARAMETERS bohr",
            );
        }
        cellUnits = coefficients.BOHR_TO_ANGSTROM;
    } else if (units === "angstrom") {
        if (alat !== null) {
            throw new Error(
                "Lattice parameters given in &SYSTEM celldm/A and CELL_PARAMETERS angstrom",
            );
        }
        cellUnits = 1.0;
    } else if (units === "alat") {
        if (alat === null) {
            throw new Error("Lattice parameters must be set in &SYSTEM for alat units");
        } else cellUnits = alat * coefficients.BOHR_TO_ANGSTROM;
    } else if (alat === null) {
        cellUnits = coefficients.BOHR_TO_ANGSTROM;
    } else {
        cellUnits = alat;
    }

    const cellParamsData = cellParamsMatch[3].trim().split("\n");

    const cellVectors = [];
    cellParamsData.forEach((line) => {
        const values = line
            .trim()
            .match(/([\d.-]+)/g)
            .map(Number)
            .map((value) => value * cellUnits);
        cellVectors.push(values);
    });

    // eslint-disable-next-line no-param-reassign
    alat = 1; // TODO: remove and change this
    return {
        cell: cellVectors,
        alat,
        units,
    };
}

/**
 * @summary Creates cell from ibrav and celldm(i) parameters
 * @param {Object} system
 * @returns {Number[][], Number, String}
 */
function ibravToCell(system) {
    let cell,
        sinab,
        sinac,
        tx,
        ty,
        tz,
        // eslint-disable-next-line no-unused-vars
        a_prime,
        u,
        v,
        v3,
        alat,
        b_over_a,
        c_over_a,
        cosab,
        cosac,
        cosbc;
    if (system.celldm[1] && system.a) {
        throw new Error("do not specify both celldm and a,b,c!");
    } else if (system.celldm[1]) {
        // celldm(x) in bohr
        // eslint-disable-next-line prefer-destructuring
        alat = system.celldm[1];
        b_over_a = system.celldm[2] || 0.0;
        c_over_a = system.celldm[3] || 0.0;
        cosab = system.celldm[4] || 0.0;
        cosac = system.celldm[5] || 0.0;
        cosbc = 0.0;
        if (system.ibrav === 14) {
            cosbc = system.celldm[4] || 0.0;
            cosac = system.celldm[5] || 0.0;
            cosab = system.celldm[6] || 0.0;
        }
    } else if (system.a) {
        // a, b, c, cosAB, cosAC, cosBC in Angstrom
        throw new Error("params_to_cell() does not yet support A/B/C/cosAB/cosAC/cosBC");
    } else {
        throw new Error("Missing celldm(1)");
    }

    switch (system.ibrav) {
        case 1:
            cell = [
                [alat, 0, 0],
                [0, alat, 0],
                [0, 0, alat],
            ];
            break;
        case 2:
            cell = [
                [-0.5 * alat, 0, 0.5 * alat],
                [0, 0.5 * alat, 0.5 * alat],
                [-0.5 * alat, 0.5 * alat, 0],
            ];
            break;
        case 3:
            cell = [
                [0.5 * alat, 0.5 * alat, 0.5 * alat],
                [-0.5 * alat, 0.5 * alat, 0.5 * alat],
                [-0.5 * alat, -0.5 * alat, 0.5 * alat],
            ];
            break;
        case -3:
            cell = [
                [-0.5 * alat, 0.5 * alat, 0.5 * alat],
                [0.5 * alat, -0.5 * alat, 0.5 * alat],
                [0.5 * alat, 0.5 * alat, -0.5 * alat],
            ];
            break;
        case 4:
            cell = [
                [alat, 0, 0],
                [-0.5 * alat, 0.5 * alat * math.sqrt(3), 0],
                [0, 0, c_over_a * alat],
            ];
            break;
        case 5:
            tx = math.sqrt((1.0 - cosab) / 2.0);
            ty = math.sqrt((1.0 - cosab) / 6.0);
            tz = math.sqrt((1 + 2 * cosab) / 3.0);
            cell = [
                [tx * alat, -ty * alat, tz * alat],
                [0, 2 * ty * alat, tz * alat],
                [-tx * alat, -ty * alat, tz * alat],
            ];
            break;
        case -5:
            u = (1.0 - cosab) / 4.0;
            v = (1.0 - 2.0 * cosab + 2.0 * math.sqrt(1 + cosab)) / 4.0;
            v3 = 2.0 * math.sqrt(v);
            cell = [
                [u * alat, v3 * alat, (0.5 - u) * alat],
                [-u * alat, v3 * alat, (0.5 + u) * alat],
                [(0.5 - v) * alat, -v3 * alat, -0.5 * alat],
            ];
            break;
        case 6:
            cell = [
                [alat, 0, 0],
                [0, b_over_a * alat, 0],
                [0, 0, c_over_a * alat],
            ];
            break;
        case 7:
            cell = [
                [alat, 0, 0],
                [0, b_over_a * alat, 0],
                [0, b_over_a * alat * cosab, math.sqrt((b_over_a * alat) ** 2 * (1 - cosab ** 2))],
            ];
            break;
        case 8:
            cell = [
                [alat, 0, 0],
                [0, b_over_a * alat, 0],
                [0, 0, c_over_a * alat],
            ];
            break;
        case 9:
            cell = [
                [0.5 * alat, (b_over_a * alat) / 2, 0],
                [-0.5 * alat, (b_over_a * alat) / 2, 0],
                [0, 0, c_over_a * alat],
            ];
            break;
        case -9:
            cell = [
                [0.5 * alat, (-b_over_a * alat) / 2, 0],
                [0.5 * alat, (b_over_a * alat) / 2, 0],
                [0, 0, c_over_a * alat],
            ];
            break;
        case 10:
            cell = [
                [alat, 0, 0],
                [alat * cosab, alat * math.sin(math.acos(cosab)), 0],
                [c_over_a * alat * cosac, 0, c_over_a * alat * math.sin(math.acos(cosac))],
            ];
            break;
        case 11:
            cell = [
                [alat * 0.5, (-alat * math.sqrt(3)) / 2, 0],
                [alat * 0.5, (alat * math.sqrt(3)) / 2, 0],
                [alat * cosab, alat * math.sqrt(1 - cosab ** 2), c_over_a * alat],
            ];
            break;
        case 12:
            sinab = math.sqrt(1.0 - cosab ** 2);
            cell = [
                [alat, 0, 0],
                [b_over_a * alat * cosab, b_over_a * alat * sinab, 0],
                [0, 0, c_over_a * alat],
            ];
            break;
        case -12:
            sinac = math.sqrt(1.0 - cosac ** 2);
            cell = [
                [alat, 0, 0],
                [0, b_over_a * alat, 0],
                [c_over_a * alat * cosac, 0, c_over_a * alat * sinac],
            ];
            break;
        case 13:
            sinab = math.sqrt(1.0 - cosab ** 2);
            cell = [
                [0.5 * alat, 0, (-c_over_a * alat) / 2],
                [b_over_a * alat * cosab, b_over_a * alat * sinab, 0],
                [0.5 * alat, 0, (c_over_a * alat) / 2],
            ];
            break;
        case 14:
            sinab = math.sqrt(1.0 - cosab ** 2);
            v3 = [
                c_over_a * cosac,
                (c_over_a * (cosbc - cosac * cosab)) / sinab,
                (c_over_a *
                    math.sqrt(
                        1 + 2 * cosbc * cosac * cosab - cosbc ** 2 - cosac ** 2 - cosab ** 2,
                    )) /
                    sinab,
            ];
            cell = [
                [alat, 0, 0],
                [b_over_a * alat * cosab, b_over_a * alat * sinab, 0],
                [v3[0] * alat, v3[1] * alat, v3[2] * alat],
            ];
            break;
        default:
            throw new Error(`ibrav = ${system.ibrav} not implemented`);
    }
    const units = ATOMIC_COORD_UNITS.cartesian;
    alat = 1; // to see if this will give the correct result
    return { cell, alat, units };
}

/**
 * @summary function to get atomic positions from the ATOMIC_POSITIONS card
 * @param {String} cards
 * @param {Number} nAtoms
 * @returns {[String], [Number], [Boolean], String}
 */
function getAtomicPositions(cards) {
    const regex =
        /ATOMIC_POSITIONS\s\(([\w]+)\)(?:\n?([!#].*)?)\n((?:\w+\s+[-?\d+.\d+\s+]*\n)*(?!\w+\s+\())?/gm;
    const match = regex.exec(cards);

    if (!match) {
        throw new Error("ATOMIC_POSITIONS section not found in input data.");
    }

    const coordSystem = match[1];
    switch (coordSystem) {
        case "angstrom":
            break;
        case "bohr":
            break;
        case "crystal":
            break;
        case "alat":
            break;
        default:
            throw new Error(`unknown coordinate system: ${coordSystem}`);
    }
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

    const units = coordSystem;

    return {
        elements,
        coordinates,
        constraints,
        units,
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
    const { ibrav } = data.system;
    const { cards } = data;

    let cell = null;
    if (ibrav === 0) {
        cell = getCellParameters(cards);
    } else {
        cell = ibravToCell(data.system);
    }
    const { elements, coordinates, constraints, units } = getAtomicPositions(cards);

    const lattice = Lattice.fromVectors({
        a: cell.cell[0],
        b: cell.cell[1],
        c: cell.cell[2],
        alat: cell.alat,
        units: cell.units,
    });
    // TODO: add lattice.type

    const basis = new ConstrainedBasis({
        elements,
        coordinates,
        units,
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
