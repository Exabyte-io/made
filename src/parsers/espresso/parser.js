import { ATOMIC_COORD_UNITS, coefficients } from "@exabyte-io/code.js/dist/constants";
import { mix } from "mixwith";

import { primitiveCell } from "../../cell/primitive_cell";
import { Lattice } from "../../lattice/lattice";
import math from "../../math";
import { MaterialParser } from "../structure";
import { FortranParserMixin } from "../utils/fortran";
import { IBRAV_TO_LATTICE_TYPE_MAP, regex } from "./settings";

export class ESPRESSOMaterialParser extends mix(MaterialParser).with(FortranParserMixin) {
    parse(content) {
        this.data = this.parseNamelists(content);
        return this.parseMaterial();
    }

    getCell() {
        return this.getCellConfig(this.data.cards);
    }

    getElements() {
        const text = this.data.cards;
        const atomicPositionsMatches = Array.from(text.matchAll(regex.atomicPositions));
        const elements = atomicPositionsMatches.map((match, index) => ({
            id: index,
            value: match[1],
        }));
        return elements;
    }

    getCoordinates() {
        const text = this.data.cards;
        const atomicPositionsMatches = Array.from(text.matchAll(regex.atomicPositions));
        const units = text.match(regex.atomicPositionsUnits)[1];
        const { scalingFactor } = this.getCoordinatesUnitsScalingFactor(units);
        const coordinates = atomicPositionsMatches.map((match, index) => ({
            id: index,
            value: [
                parseFloat(match[2]) * scalingFactor,
                parseFloat(match[3]) * scalingFactor,
                parseFloat(match[4]) * scalingFactor,
            ],
        }));
        return coordinates;
    }

    getConstraints() {
        const text = this.data.cards;
        const atomicPositionsMatches = Array.from(text.matchAll(regex.atomicPositions));
        const constraints = atomicPositionsMatches.reduce((acc, match, index) => {
            if (match[5] && match[6] && match[7]) {
                // Check if all three constraints exist
                acc.push({
                    id: index,
                    value: [match[5] === "1", match[6] === "1", match[7] === "1"],
                });
            }
            return acc;
        }, []);
        return constraints;
    }

    getUnits() {
        const text = this.data.cards;
        const units = text.match(regex.atomicPositionsUnits)[1];
        return this.getCoordinatesUnitsScalingFactor(units).units;
    }

    getName() {
        return this.data.control.title;
    }

    /**
     * @summary Return unit cell parameters from CELL_PARAMETERS card
     * @param {String} text - cards data
     * @return {{cell: Number[][], units: String}}
     */
    getCellConfig(text) {
        let cell = {};
        if (this.data.system === undefined)
            throw new Error("No &SYSTEM section found in input this.data.");
        if (this.data.system.ibrav === undefined) throw new Error("ibrav is required in &SYSTEM.");

        if (this.data.system.ibrav === 0) {
            const match = regex.cellParameters.exec(text);
            if (match) {
                const units = match[1];
                const values = match.slice(2, 11);
                // creating matrix 3 by 3 of numbers from 9 strings
                const vectors = Array.from({ length: 3 }, (_, i) =>
                    values.slice(i * 3, i * 3 + 3).map(Number),
                );
                cell = { cell: vectors, units };
                cell.type = "TRI"; // TODO: implement type detection, now defaults to TRI
                return cell;
            }
        } else {
            cell = this.ibravToCellConfig(this.data.system);
            return cell;
        }
        throw new Error("Couldn't read cell parameters");
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
     * @returns {{cell: Number[][], type: String}}
     */
    ibravToCellConfig(system) {
        const { ibrav, celldm, cosab, cosac, cosbc } = system;
        let { a, b, c } = system;
        if (celldm && a) {
            throw new Error("Both celldm and A are given");
        } else if (!celldm && !a) {
            throw new Error("Missing celldm(1)");
        }

        const type = this.ibravToCellType(ibrav);
        [a, b, c] = celldm ? this.getLatticeConstants(celldm) : [a, b, c];
        const [alpha, beta, gamma] = this.getLatticeAngles(celldm, cosbc, cosac, cosab);

        const lattice = new Lattice({
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
            type,
        });
        const cell = primitiveCell(lattice);
        return { cell, type };
    }

    /**
     * @summary Converts ibrav value to cell type according to Quantum ESPRESSO docs
     * https://www.quantum-espresso.org/Doc/INPUT_PW.html#ibrav
     * @param {Number} ibrav - ibrav parameter
     * @returns {String}
     */
    ibravToCellType(ibrav) {
        const type = IBRAV_TO_LATTICE_TYPE_MAP[ibrav];
        if (type === undefined) {
            throw new Error(`Invalid ibrav value: ${ibrav}`);
        }
        return type;
    }

    /**
     * @summary Calculates cell parameters from celldm(i) or A, B, C parameters depending on which are present. Specific to Quantum ESPRESSO.
     * @param {Number[]} celldm - celldm(i) parameters
     * @returns {Number[]}
     */
    getLatticeConstants(celldm) {
        // celldm indices shifted -1 from fortran list representation. In QE input file celldm(1) list starts with 1, but parsed starting with 0.
        const a = celldm[0] * coefficients.BOHR_TO_ANGSTROM; // celldm(1) is a in bohr
        const b = celldm[1] * celldm[0] * coefficients.BOHR_TO_ANGSTROM; // celldm(2) is b/a
        const c = celldm[2] * celldm[0] * coefficients.BOHR_TO_ANGSTROM; // celldm(3) is c/a
        return [a, b, c];
    }

    /**
     * @summary Calculates cell angles from celldm(i) or cosAB, cosAC, cosBC parameters. Specific to Quantum ESPRESSO.
     * @param {Number[]} [celldm] - celldm(i) parameters
     * @param {Number} [cosbc] - cosBC parameter
     * @param {Number} [cosac] - cosAC parameter
     * @param {Number} [cosab]   - cosAB parameter
     * @returns {Array<Number | undefined>}
     */
    getLatticeAngles(celldm, cosbc, cosac, cosab) {
        let alpha, beta, gamma;
        if (cosbc) alpha = math.acos(cosbc);
        if (cosac) beta = math.acos(cosac);
        if (cosab) gamma = math.acos(cosab);

        // Case for some of the cell types in QE docs
        // celldm indices shifted -1 from fortran list representation. In QE input file celdm(1) array starts with 1, but parsed starting with 0.
        if (celldm && celldm[3]) {
            gamma = math.acos(celldm[3]);
        }

        // Specific case for hexagonal cell in QE docs
        // celldm indices shifted -1 from fortran list representation. In QE input file celdm(1) array starts with 1, but parsed starting with 0.
        if (celldm && celldm[3] && celldm[4] && celldm[5]) {
            alpha = math.acos(celldm[3]);
            beta = math.acos(celldm[4]);
            gamma = math.acos(celldm[5]);
        }

        // Convert radians to degrees which are used in lattice definitions
        [alpha, beta, gamma] = [alpha, beta, gamma].map((x) =>
            x === undefined ? x : (x * 180) / math.PI,
        );
        return [alpha, beta, gamma];
    }

    /**
     * @summary Return units and scaling factor according to Quantum ESPRESSO docs
     * @param {String} units - units from ATOMIC_POSITIONS card
     * @returns {{units: String, scalingFactor: Number}}
     */
    getCoordinatesUnitsScalingFactor(units) {
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
        return { units: _units, scalingFactor };
    }
}
