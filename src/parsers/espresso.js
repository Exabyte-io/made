import _ from "underscore";
import s from "underscore.string";

import { ConstrainedBasis } from "../basis/constrained_basis";
// eslint-disable-next-line no-unused-vars
import { ATOMIC_COORD_UNITS, coefficients } from "../constants";
import { Lattice } from "../lattice/lattice";
import { LATTICE_TYPE } from "../lattice/types";
import math from "../math";
import { parseFortranFile } from "./fortran_parser";
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
    const regex = /&CONTROL|&SYSTEM|ATOMIC_POSITIONS/i;
    return regex.test(fileContent);
}

// /**
//  * @summary Function to convert Fortran namelist to JavaScript object.
//  * @param {String} data - Fortran string.
//  * @return {Object}
//  */
// function parseFortranNameList(data) {
//     const blocks = data.split(/\/\n+/);
//     const result = {};
//
//     blocks.forEach((block) => {
//         const lines = block.split("\n");
//         const blockName = lines[0].trim().toLowerCase();
//
//         if (blockName.startsWith("&")) {
//             const blockContent = {};
//             for (let i = 1; i < lines.length; i++) {
//                 const line = lines[i].trim();
//                 if (line) {
//                     // Check if the line contains celldm(i)
//                     const celldmMatch = line.match(/celldm\((\d+)\)\s*=\s*(\d.\d+).*(#.*)?/);
//                     if (celldmMatch) {
//                         const index = Number(celldmMatch[1]);
//                         const value = Number(celldmMatch[2]);
//                         if (!blockContent.celldm) {
//                             blockContent.celldm = {}; // Initialize if not present
//                         }
//                         blockContent.celldm[index] = value; // Store value at index
//                     } else {
//                         const [key, value] = line.split("=").map((str) => str.trim());
//                         // convert Fortran boolean to JavaScript boolean
//                         if (value === ".true.") {
//                             blockContent[key] = true;
//                         } else if (value === ".false.") {
//                             blockContent[key] = false;
//                         } else if (Number.isNaN(Number(value))) {
//                             // if value is not a number, keep it as a string
//                             blockContent[key] = value.replace(/'/g, ""); // remove single quotes
//                         } else {
//                             blockContent[key] = Number(value); // convert to number
//                         }
//                     }
//                 }
//             }
//             result[blockName.slice(1)] = blockContent;
//         } else {
//             result.cards = block;
//         }
//     });
//
//     return result;
// }

/**
 * @summary Read unit cell parameters from CELL_PARAMETERS card
 * @param {Object} cell
 * @param {Number} alat
 * @returns {cell: Number[][], alat: Number, units: String}
 */
function getCellParameters(cell, alat = null) {
    let { units } = cell;
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

    Object.keys(cell.cell).forEach((key) => {
        if (Array.isArray(cell.cell[key])) {
            cell.cell[key] = cell.cell[key].map((value) => value * cellUnits);
        }
    });

    // eslint-disable-next-line no-param-reassign
    alat = 1; // Supplied data already taken in consideration so here it set to unity
    // eslint-disable-next-line no-param-reassign
    units = ATOMIC_COORD_UNITS.angstrom; // Supplied data already taken in consideration
    cell.alat = alat;
    return cell;
}

// TODO: move lattice logic to Lattice
/**
 * @summary Creates cell from ibrav and celldm(i) parameters
 * // https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm224
 * @param {Object} system ???
 * @returns {Number[][], Number, String, String}
 */
function ibravToCell(system) {
    let alat = 1.0;
    let cell,
        type,
        sinab,
        sinac,
        tx,
        ty,
        tz,
        a_prime,
        u,
        v,
        v3,
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
        alat = system.celldm[1] * coefficients.BOHR_TO_ANGSTROM; // this value given in Bohr (a.u.) as per QE specs
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
            type = LATTICE_TYPE.CUB;
            break;
        case 2:
            cell = [
                [-0.5 * alat, 0, 0.5 * alat],
                [0, 0.5 * alat, 0.5 * alat],
                [-0.5 * alat, 0.5 * alat, 0],
            ];
            type = LATTICE_TYPE.FCC;
            break;
        case 3:
            cell = [
                [0.5 * alat, 0.5 * alat, 0.5 * alat],
                [-0.5 * alat, 0.5 * alat, 0.5 * alat],
                [-0.5 * alat, -0.5 * alat, 0.5 * alat],
            ];
            type = LATTICE_TYPE.BCC;
            break;
        case -3:
            cell = [
                [-0.5 * alat, 0.5 * alat, 0.5 * alat],
                [0.5 * alat, -0.5 * alat, 0.5 * alat],
                [0.5 * alat, 0.5 * alat, -0.5 * alat],
            ];
            type = LATTICE_TYPE.BCC;
            break;
        case 4:
            cell = [
                [alat, 0, 0],
                [-0.5 * alat, 0.5 * alat * math.sqrt(3), 0],
                [0, 0, c_over_a * alat],
            ];
            type = LATTICE_TYPE.HEX;
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
            type = LATTICE_TYPE.TRI; // FIXME: NOT TRUE!
            break;
        case -5:
            ty = math.sqrt((1.0 - cosab) / 6.0);
            tz = math.sqrt((1 + 2 * cosab) / 3.0);
            a_prime = alat / math.sqrt(3);
            u = tz - 2 * math.sqrt(2) * ty;
            v = tz + math.sqrt(2) * ty;
            cell = [
                [u * a_prime, v * a_prime, v * a_prime],
                [v * a_prime, u * a_prime, v * a_prime],
                [v * a_prime, v * a_prime, u * a_prime],
            ];
            type = LATTICE_TYPE.TRI; // FIXME: NOT TRUE!
            break;
        case 6:
            cell = [
                [alat, 0, 0],
                [0, alat, 0],
                [0, 0, c_over_a * alat],
            ];
            type = LATTICE_TYPE.TET;
            break;
        case 7:
            cell = [
                [0.5 * alat, -0.5 * alat, c_over_a],
                [0.5 * alat, 0.5 * alat, c_over_a],
                [-0.5 * alat, -0.5 * alat, c_over_a * 0.5 * alat],
            ];
            type = LATTICE_TYPE.TET;
            break;
        case 8:
            cell = [
                [alat, 0, 0],
                [0, b_over_a * alat, 0],
                [0, 0, c_over_a * alat],
            ];
            type = LATTICE_TYPE.ORC;
            break;
        case 9:
            cell = [
                [0.5 * alat, (b_over_a * alat) / 2, 0],
                [-0.5 * alat, (b_over_a * alat) / 2, 0],
                [0, 0, c_over_a * alat],
            ];
            type = LATTICE_TYPE.ORCC;
            break;
        case -9:
            cell = [
                [0.5 * alat, (-b_over_a * alat) / 2, 0],
                [0.5 * alat, (b_over_a * alat) / 2, 0],
                [0, 0, c_over_a * alat],
            ];
            type = LATTICE_TYPE.ORCC;
            break;
        case 10:
            cell = [
                [0.5 * alat, 0, 0.5 * c_over_a],
                [0.5 * alat, 0.5 * b_over_a * alat, 0],
                [0, 0.5 * b_over_a * alat, 0.5 * c_over_a * alat],
            ];
            type = LATTICE_TYPE.ORCF;
            break;
        case 11:
            cell = [
                [0.5 * alat, 0.5 * b_over_a * alat, 0.5 * c_over_a * alat],
                [-0.5 * alat, 0.5 * b_over_a * alat, 0.5 * c_over_a * alat],
                [-0.5 * alat, -0.5 * b_over_a * alat, 0.5 * c_over_a * alat],
            ];
            type = LATTICE_TYPE.ORCI;
            break;
        case 12:
            sinab = math.sqrt(1.0 - cosab ** 2);
            cell = [
                [alat, 0, 0],
                [b_over_a * alat * cosab, b_over_a * alat * sinab, 0],
                [0, 0, c_over_a * alat],
            ];
            type = LATTICE_TYPE.MCL; // Monoclinic P, unique axis c
            break;
        case -12:
            sinac = math.sqrt(1.0 - cosac ** 2);
            cell = [
                [alat, 0, 0],
                [0, b_over_a * alat, 0],
                [c_over_a * alat * cosac, 0, c_over_a * alat * sinac],
            ];
            type = LATTICE_TYPE.MCL; // Monoclinic P, unique axis b
            break;
        case 13:
            sinab = math.sqrt(1.0 - cosab ** 2);
            cell = [
                [0.5 * alat, 0, -0.5 * c_over_a * alat],
                [b_over_a * alat * cosab, b_over_a * alat * sinab, 0],
                [0.5 * alat, 0, 0.5 * c_over_a * alat],
            ];
            type = LATTICE_TYPE.MCLC; // Monoclinic P, base centered
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
            type = LATTICE_TYPE.TRI;
            break;
        default:
            throw new Error(`ibrav = ${system.ibrav} not implemented`);
    }
    const units = ATOMIC_COORD_UNITS.angstrom; // This is passed to the
    alat = 1.0;
    return { cell, alat, units, type };
}

/**
 * Parses QE .in file into a Material config object.
 * @param {String} fileContent - contents of the .in file
 * @return {Object} Material config.
 */
function fromEspressoFormat(fileContent) {
    const data = parseFortranFile(fileContent);
    // const nSpecies = data.system.ntyp;
    // const nAtoms = data.system.nat;
    const { ibrav } = data.system;
    let cell = null;

    // Check for errors
    if (data.system === undefined) {
        throw new Error("No SYSTEM section found in input data.");
    }
    if (ibrav === undefined) {
        throw new Error("ibrav is required in &SYSTEM.");
    }

    if (ibrav === 0) {
        cell = getCellParameters(data.cell);
    } else {
        // Create unit cell from given ibrav (and celldm(i) or A) with algorithm
        cell = ibravToCell(data.system);
    }
    const { elements, coordinates, constraints, units } = data;
    const lattice = Lattice.fromVectors({
        a: cell.cell[0],
        b: cell.cell[1],
        c: cell.cell[2],
        alat: cell.alat,
        units: cell.units,
    });
    if (cell.type) lattice.type = cell.type;

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
