import { coefficients } from "@exabyte-io/code.js/dist/constants";
import _ from "underscore";
import s from "underscore.string";

import { ConstrainedBasis } from "../basis/constrained_basis";
import { Lattice } from "../lattice/lattice";
import { LATTICE_TYPE } from "../lattice/types";
import math from "../math";
import { parseFortranFile } from "./fortran_parser";
import { regex } from "./utils";
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
    const espressoRegex = /&CONTROL|&SYSTEM|ATOMIC_POSITIONS/i;
    return espressoRegex.test(fileContent);
}

/**
 * @summary Read unit cell parameters from CELL_PARAMETERS card
 * @param {String} text - cards data
 * @param {Number} alat - lattice parameter
 * @returns {{vectors: number[][], alat: Number, units: string}}
 */
function getCellParameters(text, alat = null) {
    let vectors = [];
    let units;

    const cellParameters = text.match(regex.cellParameters);

    if (cellParameters) {
        vectors = [
            [
                parseFloat(cellParameters[2]),
                parseFloat(cellParameters[3]),
                parseFloat(cellParameters[4]),
            ],
            [
                parseFloat(cellParameters[5]),
                parseFloat(cellParameters[6]),
                parseFloat(cellParameters[7]),
            ],
            [
                parseFloat(cellParameters[8]),
                parseFloat(cellParameters[9]),
                parseFloat(cellParameters[10]),
            ],
        ];
        // eslint-disable-next-line prefer-destructuring
        units = cellParameters[1];
    }

    if (units === "bohr") {
        if (alat !== null) {
            throw new Error(
                "Lattice parameters given in &SYSTEM celldm/A and CELL_PARAMETERS bohr",
            );
        }
    } else if (units === "angstrom") {
        if (alat !== null) {
            throw new Error(
                "Lattice parameters given in &SYSTEM celldm/A and CELL_PARAMETERS angstrom",
            );
        }
    } else if (units === "alat") {
        if (alat === null) {
            throw new Error("Lattice parameters must be set in &SYSTEM for alat units");
        }
    }

    // eslint-disable-next-line no-param-reassign
    alat = 1; // Supplied data already taken in consideration so here it set to unity
    return { vectors, alat, units };
}

/**
 * @summary Read atomic positions from ATOMIC_POSITIONS card
 * @param {String} text - cards data
 * @returns {{elements: Object[], coordinates: Object[], constraints: Object[], units: String}}
 */
function getAtomicPositions(text) {
    const atomicSpeciesMatches = Array.from(text.matchAll(regex.atomicSpecies));
    // eslint-disable-next-line no-unused-vars
    const atomicSpecies = atomicSpeciesMatches.map((match) => ({
        element: match[1],
        mass: parseFloat(match[2]),
        potential: match[3],
    }));
    const atomicPositionsMatches = Array.from(text.matchAll(regex.atomicPositions));
    const elements = atomicPositionsMatches.map((match, index) => ({
        id: index,
        value: match[1],
    }));
    const coordinates = atomicPositionsMatches.map((match, index) => ({
        id: index,
        value: [parseFloat(match[2]), parseFloat(match[3]), parseFloat(match[4])],
    }));
    const constraints = atomicPositionsMatches
        .filter((match) => match[5] && match[6] && match[7]) // Check if all three constrtaints exist
        .map((match, index) => ({
            id: index,
            value: [match[5] === "1", match[6] === "1", match[7] === "1"],
        }));
    // eslint-disable-next-line prefer-destructuring
    const units = text.match(regex.atomicPositionsUnits)[1];

    return { elements, coordinates, constraints, units };
}

/**
 * QE docs: https://www.quantum-espresso.org/Doc/INPUT_PW.html#ibrav
 * "If ibrav /= 0, specify EITHER [ celldm(1)-celldm(6) ]
 *   OR [ A, B, C, cosAB, cosAC, cosBC ]
 *   but NOT both. The lattice parameter "alat" is set to
 *   alat = celldm(1) (in a.u.) or alat = A (in Angstrom);"
 *
 * @summary Creates cell from ibrav and celldm(i) parameters
 * @param {Object} system - The system parameters from &SYSTEM namelist
 * @param {Number} system.ibrav - ibrav parameter
 * @param {Number[]} system.celldm - celldm parameters, 0 - a (bohr), 1 - b/a (unitless), 2 - c/a (unitless), 3,4,5 - cosines
 * @param {Number} [system.a] - a parameter in angstroms
 * @param {Number} [system.b] - b parameter in angstroms
 * @param {Number} [system.c] - c parameter in angstroms
 * @param {Number} [system.cosab] - cosab parameter
 * @param {Number} [system.cosac] - cosac parameter
 * @param {Number} [system.cosbc] - cosbc parameter
 * @returns {{vectors: number[][], alat: Number, units: string}}
 */
function ibravToCell(system) {
    const { ibrav, celldm, a, b, c, cosab, cosac, cosbc } = system;
    const alat = 1;
    let type, alpha, beta, gamma;
    let _a = celldm[0] ? celldm[0] : a;
    let _b = celldm[1] ? celldm[1] * celldm[0] : b; // celldm[1] is b/a
    let _c = celldm[2] ? celldm[2] * celldm[0] : c; // celldm[2] is c/a

    // celdm(1) = celldm[0]: index begins from 0
    if (celldm[0] && a) {
        throw new Error("Both celldm and A are given");
    } else if (celldm[0]) {
        // celldm(x) in bohr - from QE docs
        _a *= coefficients.BOHR_TO_ANGSTROM;
        _b *= coefficients.BOHR_TO_ANGSTROM;
        _c *= coefficients.BOHR_TO_ANGSTROM;
    } else if (a) {
        // a, b, c, cosAB, cosAC, cosBC in Angstrom - from QE docs
        alpha = math.acos(cosbc);
        beta = math.acos(cosac);
        gamma = math.acos(cosab);
    } else {
        throw new Error("Missing celldm(1)");
    }

    switch (ibrav) {
        case 1:
            type = LATTICE_TYPE.CUB;
            break;
        case 2:
            type = LATTICE_TYPE.FCC;
            break;
        case 3:
            type = LATTICE_TYPE.BCC;
            break;
        case -3:
            type = LATTICE_TYPE.BCC;
            break;
        case 4:
            type = LATTICE_TYPE.HEX;
            break;
        case 5:
            type = LATTICE_TYPE.RHL;
            break;
        case -5:
            type = LATTICE_TYPE.RHL;
            break;
        case 6:
            type = LATTICE_TYPE.TET;
            break;
        case 7:
            type = LATTICE_TYPE.TET;
            break;
        case 8:
            type = LATTICE_TYPE.ORC;
            break;
        case 9:
            type = LATTICE_TYPE.ORCC;
            break;
        case -9:
            type = LATTICE_TYPE.ORCC;
            break;
        case 10:
            type = LATTICE_TYPE.ORCF;
            break;
        case 11:
            type = LATTICE_TYPE.ORCI;
            break;
        case 12:
            type = LATTICE_TYPE.MCL; // Monoclinic P, unique axis c
            break;
        case -12:
            type = LATTICE_TYPE.MCL; // Monoclinic P, unique axis b
            break;
        case 13:
            type = LATTICE_TYPE.MCLC; // Monoclinic P, base centered
            break;
        case 14:
            type = LATTICE_TYPE.TRI;
            break;
        default:
            throw new Error(`ibrav = ${system.ibrav} not implemented`);
    }

    const { vectors } = Lattice.vectorsFromType(type, _a, _b, _c, alpha, beta, gamma);
    return { vectors, alat, units: "angstrom", type };
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
    let cell = {};
    let lattice = {};

    // Check for errors
    if (data.system === undefined) {
        throw new Error("No SYSTEM section found in input data.");
    }
    if (ibrav === undefined) {
        throw new Error("ibrav is required in &SYSTEM.");
    }

    if (ibrav === 0) {
        // Create unit cell from given cell parameters in CELL_PARAMETERS card
        cell = getCellParameters(data.cards);
    } else {
        // Create unit cell from given ibrav (and celldm(i) or A) with algorithm
        cell = ibravToCell(data.system);
    }

    lattice = Lattice.fromVectors({
        a: cell.vectors[0],
        b: cell.vectors[1],
        c: cell.vectors[2],
        alat: cell.alat,
        units: cell.units,
        type: cell.type,
    });

    const { elements, coordinates, constraints, units } = getAtomicPositions(data.cards);

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
