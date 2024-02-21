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
        this.data = this.fortranParseNamelists(content);
        return this.parseMaterial();
    }

    /**
     * @summary Return unit cell parameters from CELL_PARAMETERS card
     * @returns {{cell: Number[][], units: String}}
     */
    getCell() {
        const text = this.data.cards;
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
                // TODO: implement type detection, now defaults to TRI
                cell.type = "TRI";
                return cell;
            }
        } else {
            cell = this.ibravToCellConfig();
            return cell;
        }
        throw new Error("Couldn't read cell parameters");
    }

    /**
     * @summary Return elements from ATOMIC_SPECIES card
     * @returns {{id: Number, value: String}[]}
     */
    getElements() {
        const text = this.data.cards;
        const atomicPositionsMatches = Array.from(text.matchAll(regex.atomicPositions));
        return atomicPositionsMatches.map((match, index) => ({
            id: index,
            value: match[1],
        }));
    }

    /**
     * @summary Return atomic positions from ATOMIC_POSITIONS card
     * @returns {{id: Number, value: Number[]}[]}
     */
    getCoordinates() {
        const text = this.data.cards;
        const atomicPositionsMatches = Array.from(text.matchAll(regex.atomicPositions));
        const { scalingFactor } = this.getCoordinatesUnitsScalingFactor();
        return atomicPositionsMatches.map((match, index) => ({
            id: index,
            value: match.slice(2, 5).map((value) => parseFloat(value) * scalingFactor),
        }));
    }

    /**
     * @summary Return atomic constraints from ATOMIC_POSITIONS card
     * @returns {{id: Number, value: Boolean[]}[]}
     */
    getConstraints() {
        const text = this.data.cards;
        const atomicPositionsMatches = Array.from(text.matchAll(regex.atomicPositions));

        const constraints = atomicPositionsMatches.reduce((acc, match, index) => {
            const value = match
                .slice(5, 8)
                .filter((constraint) => constraint !== undefined)
                .map((constraint) => constraint === "1"); // expect only 0 or 1 as valid values

            acc.push({
                id: index,
                value,
            });

            return acc;
        }, []);

        // If all constraints are empty, return an empty array
        if (constraints.every((constraint) => constraint.value.length === 0)) {
            return [];
        }

        return constraints;
    }

    /**
     * @summary Return atomic coordinates units from ATOMIC_POSITIONS card
     * @returns {String}
     */
    getUnits() {
        return this.getCoordinatesUnitsScalingFactor().units;
    }

    /**
     * @summary Return material name from CONTROL card
     * If not present, later will be generated from the formula in materialConfig object
     * @returns {String}
     */
    getName() {
        return this.data.control.title;
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
     * @returns {{cell: Number[][], type: String}}
     */
    ibravToCellConfig() {
        const { system } = this.data;
        const { celldm } = system;
        let { a, b, c } = system;
        if (celldm && a) {
            throw new Error("Both celldm and A are given");
        } else if (!celldm && !a) {
            throw new Error("Missing celldm(1)");
        }

        const type = this.ibravToCellType();
        [a, b, c] = celldm ? this.getLatticeConstants() : [a, b, c];
        const [alpha, beta, gamma] = this.getLatticeAngles();

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
     * @returns {String}
     */
    ibravToCellType() {
        const { ibrav } = this.data.system;
        const type = IBRAV_TO_LATTICE_TYPE_MAP[ibrav];
        if (type === undefined) {
            throw new Error(`Invalid ibrav value: ${ibrav}`);
        }
        return type;
    }

    /**
     * @summary Calculates cell parameters from celldm(i) or A, B, C parameters depending on which are present. Specific to Quantum ESPRESSO.
     * @returns {Number[]}
     * */
    getLatticeConstants() {
        const { celldm } = this.data.system;
        // celldm indices shifted -1 from fortran list representation. In QE input file celldm(1) list starts with 1, but parsed starting with 0.
        const a = celldm[0] * coefficients.BOHR_TO_ANGSTROM; // celldm(1) is a in bohr
        const b = celldm[1] * celldm[0] * coefficients.BOHR_TO_ANGSTROM; // celldm(2) is b/a
        const c = celldm[2] * celldm[0] * coefficients.BOHR_TO_ANGSTROM; // celldm(3) is c/a
        return [a, b, c];
    }

    /**
     * @summary Calculates cell angles from celldm(i) or cosAB, cosAC, cosBC parameters. Specific to Quantum ESPRESSO.
     * @returns {Array<Number | undefined>}
     * */
    getLatticeAngles() {
        const { celldm, cosbc, cosac, cosab } = this.data.system;
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
     * @summary Return units and scaling factor according to Quantum ESPRESSO 7.2 docs
     * @returns {{units: String, scalingFactor: Number}}
     */
    getCoordinatesUnitsScalingFactor() {
        const units = this.data.cards.match(regex.atomicPositionsUnits)[1];
        let scalingFactor = 1.0;
        let _units;
        switch (units) {
            case "alat":
                _units = ATOMIC_COORD_UNITS.crystal;
                break;
            case "bohr":
                scalingFactor = coefficients.BOHR_TO_ANGSTROM;
                _units = ATOMIC_COORD_UNITS.cartesian;
                break;
            case "angstrom":
                _units = ATOMIC_COORD_UNITS.cartesian;
                break;
            case "crystal":
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
