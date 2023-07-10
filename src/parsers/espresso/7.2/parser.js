import { ATOMIC_COORD_UNITS, coefficients } from "@exabyte-io/code.js/dist/constants";

import { primitiveCell } from "../../../cell/primitive_cell";
import { Lattice } from "../../../lattice/lattice";
import math from "../../../math";
import { FortranParser } from "../../../utils/parsers/fortran";
import { IBRAV_TO_LATTICE_TYPE_MAP, regex } from "./settings";

export class EspressoParser {
    constructor(content) {
        this.content = content;
        this.intermediateFormat = this.getIntermediateFormat();
    }

    getIntermediateFormat() {
        const data = FortranParser.parse(this.content); // using static method
        const intermediateData = data;

        intermediateData.atomicSpecies = Array.from(data.cards.matchAll(regex.atomicSpecies)).map(
            (match) => match.slice(1, match.length),
        ); // This will only take the captured groups

        intermediateData.atomicPositions = Array.from(
            data.cards.matchAll(regex.atomicPositions),
        ).map((match) => match.slice(1, match.length));

        intermediateData.cellParameters = Array.from(data.cards.matchAll(regex.cellParameters)).map(
            (match) => match.slice(1, match.length),
        );

        delete intermediateData.cards;

        console.log("intermediateFormat:", intermediateData);
        return intermediateData;
    }

    serialize() {
        const { atomicSpecies } = this.intermediateFormat();
        const { atomicPositions } = this.intermediateFormat();
        const { cellParameters } = this.intermediateFormat();
        if (!cellParameters) {
            throw new Error("Couldn't read cell parameters");
            // TODO: add cell generation from ibrav etc...
        }
        return { atomicSpecies, atomicPositions, cellParameters };
    }
}

/**
 * @summary checks if the given fileContent is in the format of a Quantum ESPRESSO .in file.
 * @param {String} fileContent
 * @return {Boolean}
 */
function isEspressoFormat(fileContent) {
    const espressoRegex = regex.espressoFingerprint;
    return espressoRegex.test(fileContent);
}

/**
 * @summary Return unit cell parameters from CELL_PARAMETERS card
 * @param {String} text - cards data
 * @return {{vectors: Number[][], units: String}}
 */
function getCellConfig(text) {
    const match = regex.cellParameters.exec(text);
    if (match) {
        const { units } = match[1];
        const values = match.slice(2, 11);
        const vectors = Array.from({ length: 3 }, (_, i) =>
            values.slice(i * 3, i * 3 + 3).map(Number),
        );
        return { vectors, units };
    }
    throw new Error("Couldn't read cell parameters");
}

/**
 * @summary Return units and scaling factor according to Quantum Espresso docs
 * @param {String} units - units from ATOMIC_POSITIONS card
 * @returns {{_units: String, scalingFactor: Number}}
 */
function convertFromEspressoUnits(units) {
    let _units, scalingFactor;
    switch (units) {
        case "alat":
            scalingFactor = 1.0;
            _units = ATOMIC_COORD_UNITS.crystal;
            break;
        case "bohr":
            scalingFactor = coefficients.BOHR_TO_ANGSTROM;
            _units = ATOMIC_COORD_UNITS.cartesian;
            break;
        case "angstrom":
            scalingFactor = 1.0;
            _units = ATOMIC_COORD_UNITS.cartesian;
            break;
        case "crystal":
            scalingFactor = 1.0;
            _units = ATOMIC_COORD_UNITS.crystal;
            break;
        case "crystal_sg":
            throw new Error("crystal_sg not supported yet");
        default:
            throw new Error(`Units ${units} not supported`);
    }
    return { _units, scalingFactor };
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
    const units = text.match(regex.atomicPositionsUnits)[1];
    const { _units, scalingFactor } = convertFromEspressoUnits(units);

    const elements = atomicPositionsMatches.map((match, index) => ({
        id: index,
        value: match[1],
    }));
    const coordinates = atomicPositionsMatches.map((match, index) => ({
        id: index,
        value: [
            parseFloat(match[2]) * scalingFactor,
            parseFloat(match[3]) * scalingFactor,
            parseFloat(match[4]) * scalingFactor,
        ],
    }));
    const constraints = atomicPositionsMatches
        .filter((match) => match[5] && match[6] && match[7]) // Check if all three constraints exist
        .map((match, index) => ({
            id: index,
            value: [match[5] === "1", match[6] === "1", match[7] === "1"],
        }));

    return { elements, coordinates, constraints, units: _units };
}

/**
 * @summary Converts ibrav value to cell type according to Quantum Espresso docs
 * https://www.quantum-espresso.org/Doc/INPUT_PW.html#ibrav
 * @param {Number} ibrav - ibrav parameter
 * @returns {String}
 */
function ibravToCellType(ibrav) {
    const type = IBRAV_TO_LATTICE_TYPE_MAP[ibrav];
    if (type === undefined) {
        throw new Error(`Invalid ibrav value: ${ibrav}`);
    }
    return type;
}

/**
 * @summary Calculates cell parameters from celldm(i) or A, B, C parameters depending on which are present. Specific to Quantum Espresso.
 * @param {Number[]} [celldm] - celldm(i) parameters
 * @param {Number} [a] - A parameter
 * @param {Number} [b] - B parameter
 * @param {Number} [c] - C parameter
 * @returns {Number[]}
 */
function getCellConstants(celldm, a, b, c) {
    let _a = celldm ? celldm[1] : a;
    let _b = celldm ? celldm[2] * celldm[1] : b; // celldm[2] is b/a
    let _c = celldm ? celldm[3] * celldm[1] : c; // celldm[3] is c/a
    if (celldm) {
        [_a, _b, _c] = [_a, _b, _c].map((x) => x * coefficients.BOHR_TO_ANGSTROM);
    }
    return [_a, _b, _c];
}

/**
 * @summary Calculates cell angles from celldm(i) or cosAB, cosAC, cosBC parameters. Specific to Quantum Espresso.
 * @param {Number[]} [celldm] - celldm(i) parameters
 * @param {Number} [cosbc] - cosBC parameter
 * @param {Number} [cosac] - cosAC parameter
 * @param {Number} [cosab]   - cosAB parameter
 * @returns {Number[]}
 */
function getCellAngles(celldm, cosbc, cosac, cosab) {
    let alpha, beta, gamma;
    if (cosbc) alpha = math.acos(cosbc);
    if (cosac) beta = math.acos(cosac);
    if (cosab) gamma = math.acos(cosab);

    // Case for some of the cell types in QE docs
    if (celldm && celldm[4]) {
        gamma = math.acos(celldm[4]);
    }

    // Specific case for hexagonal cell in QE docs
    if (celldm && celldm[4] && celldm[5] && celldm[6]) {
        alpha = math.acos(celldm[4]);
        beta = math.acos(celldm[5]);
        gamma = math.acos(celldm[6]);
    }

    // Convert radians to degrees which are used in lattice definitions
    [alpha, beta, gamma] = [alpha, beta, gamma].map((x) =>
        x === undefined ? x : (x * 180) / math.PI,
    );
    return [alpha, beta, gamma];
}

/**
 * @summary Returns cell config from ibrav and celldm(i) parameters
 *
 * QE docs: https://www.quantum-espresso.org/Doc/INPUT_PW.html#ibrav
 * "If ibrav /= 0, specify EITHER [ celldm(1)-celldm(6) ]
 *   OR [ A, B, C, cosAB, cosAC, cosBC ]
 *   but NOT both. The lattice parameter "alat" is set to
 *   alat = celldm(1) (in a.u.) or alat = A (in Angstrom);"
 *
 * @param {Object} system - The system parameters from &SYSTEM namelist
 * @param {Number} system.ibrav - ibrav parameter
 * @param {Number[]} [system.celldm] - celldm parameters
 * @param {Number} [system.a] - A parameter in angstroms
 * @param {Number} [system.b] - B parameter in angstroms
 * @param {Number} [system.c] - C parameter in angstroms
 * @param {Number} [system.cosab] - cosAB parameter
 * @param {Number} [system.cosac] - cosAC parameter
 * @param {Number} [system.cosbc] - cosBC parameter
 * @returns {{vectors: Number[][], type: String}}
 */
function ibravToCellConfig(system) {
    const { ibrav, celldm, a, b, c, cosab, cosac, cosbc } = system;
    if (celldm && a) {
        throw new Error("Both celldm and A are given");
    } else if (!celldm && !a) {
        throw new Error("Missing celldm(1)");
    }

    if (celldm) celldm.unshift(null); // Bump indices by 1 to make it easier to understand the logic. In QE input file celdm(1) array starts with 1, but parsed starting with 0.

    const type = ibravToCellType(ibrav);
    const [_a, _b, _c] = getCellConstants(celldm, a, b, c);
    const [alpha, beta, gamma] = getCellAngles(celldm, cosbc, cosac, cosab);

    const lattice = new Lattice({
        a: _a,
        b: _b,
        c: _c,
        alpha,
        beta,
        gamma,
        type,
    });
    const vectors = primitiveCell(lattice);
    return { vectors, type };
}

/**
 * Parses QE .in file into a Material config object.
 * @param {String} fileContent - contents of the .in file
 * @return {{cell: Object, elements: Object[], coordinates: Object[], constraints: Object[], units: String, name: String}}
 */
function fromEspressoFormat(fileContent) {
    const data = EspressoParser.parse(fileContent);

    if (data.system === undefined) throw new Error("No &SYSTEM section found in input data.");
    if (data.system.ibrav === undefined) throw new Error("ibrav is required in &SYSTEM.");

    let cell = {};
    if (data.system.ibrav === 0) {
        cell = getCellConfig(data.cards);
        cell.type = Lattice.typeFromVectors(cell.vectors); // not implemented yet, defaults to TRI
    } else {
        cell = ibravToCellConfig(data.system);
    }

    const { elements, coordinates, constraints, units } = getAtomicPositions(data.cards);
    const name = data.control.title;
    return { cell, elements, coordinates, constraints, units, name };
}

export default { fromEspressoFormat, isEspressoFormat };
