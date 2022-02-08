import _ from "underscore";
import s from "underscore.string";

import { Basis } from "../basis/basis";
import { ConstrainedBasis } from "../basis/constrained_basis";
import { Lattice } from "../lattice/lattice";
import math from "../math";
import { InvalidLineError } from "./errors";
import { CombinatorialBasis } from "./xyz_combinatorial_basis";

// Regular expression for an XYZ line with atomic constraints, eg. Si    0.000000    0.500000    0.446678 1 1 1`
const XYZ_LINE_REGEX =
    /[A-Z][a-z]?\s+((-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)(\s+)?(\s+[0-1]\s+[0-1]\s+[0-1](\s+)?)?)$/;

/**
 * Validates XYZ file's line. Line should be in following format "Si 0.5 0.5 0.5".
 * Raises an error if line is in wrong format.
 * @param xyzLine {String}
 * @param index {Number}
 */
export function validateLine(xyzLine, index) {
    const words = xyzLine.split(" ").filter((x) => x.trim() !== "");

    if (!xyzLine.match(XYZ_LINE_REGEX)) {
        throw new InvalidLineError(index, xyzLine);
    }

    const coordinates = [parseFloat(words[1]), parseFloat(words[2]), parseFloat(words[3])];
    coordinates.forEach((num, i) => {
        if (_.isNaN(num)) {
            throw new Error(`Coordinates should be a number. Possible error in ${i} coordinate`);
        }
    });
}

/**
 * Function checks the number of atoms listed in the BasisText editor.
 * If the number of atoms > Basis.nonPeriodicMaxAtomsCount, then an error is thrown.
 *
 * @param {String} xyzTxt
 */
export function validateNumberOfAtoms(xyzTxt) {
    const xyzArray = xyzTxt.split(/\r?\n/);
    const nAtoms = xyzArray.length;
    if (nAtoms > Basis.nonPeriodicMaxAtomsCount) {
        throw new Error("Non-Periodic max atom limit exceeded");
    }
}

/**
 * Validates that passed string is well-formed XYZ file. Function returns an error code associated with a valid (0) or
 * invalid (1001) format.
 * @param xyzTxt {String}
 * @return {Number}
 */
export function validateFormat(xyzTxt) {
    try {
        s(xyzTxt)
            .trim()
            .lines()
            .filter((x) => x.trim() !== "")
            .forEach(validateLine);
    } catch (err) {
        return 1001;
    }
    return 0;
}

/**
 * Validates that passed string is a valid XYZ file based on format and number of atoms present in the structure.
 * Function returns a number associated with whether the xyz basis is valid.
 * @param xyzTxt {String}
 * @return {Number}
 */
export function validateAll(xyzTxt) {
    let error = validateFormat(xyzTxt);
    try {
        validateNumberOfAtoms(xyzTxt);
    } catch (err) {
        error = 2001;
    }
    return error;
}

/**
 * Parses XYZ line and returns an object.
 * @param {String} line - line of text
 * @return {Object}
 */
function _parseXYZLineAsWords(line) {
    const words = s.words(line);
    return {
        element: words[0],
        coordinates: [+words[1], +words[2], +words[3]],
        // Below maps zero values to false (atom is fixed) and non-zero values to true (atom is moving)
        constraints: [+words[4], +words[5], +words[6]].map((e) => parseInt(e, 10) !== 0),
    };
}

/**
 * Parse XYZ text for basis.
 * @param txt {String} Text
 * @param units {String} Coordinate units
 * @param cell {Array} Basis Cell
 * @return {Object}
 */
function toBasisConfig(txt, units = "angstrom", cell = Basis.defaultCell) {
    const lines = s(txt).trim().lines();
    const listOfObjects = _.map(lines, _parseXYZLineAsWords);

    return {
        // using concat below to avoid modifying listOfObjects in place with map
        elements: [].concat(listOfObjects).map((elm, idx) => {
            return {
                id: idx,
                value: elm.element,
            };
        }),
        coordinates: [].concat(listOfObjects).map((elm, idx) => {
            return {
                id: idx,
                value: elm.coordinates,
            };
        }),
        units,
        cell,
        constraints: [].concat(listOfObjects).map((elm, idx) => {
            return {
                id: idx,
                value: elm.constraints,
            };
        }),
    };
}

/**
 * Create XYZ from Basis class instance.
 * @param basisClsInstance {ConstrainedBasis} Basis class instance.
 * @param printFormat {String} Output format for coordinates.
 * @param skipRounding {Boolean} Whether to round the numbers (ie. to avoid negative zeros).
 * @return {String} Basis string in XYZ format
 */
function fromBasis(basisClsInstance, printFormat = "%9.5f", skipRounding = false) {
    const clsInstance = basisClsInstance;
    const XYZArray = [];
    clsInstance._elements.array.forEach((item, idx) => {
        // assume that _elements and _coordinates are indexed equivalently
        const element = s.sprintf("%-3s", item);
        const coordinates = clsInstance
            .getCoordinateByIndex(idx)
            .map((x) =>
                s.sprintf(printFormat, skipRounding ? x : math.precise(math.roundToZero(x))),
            );
        const constraints = clsInstance.constraints
            ? clsInstance.AtomicConstraints.getAsStringByIndex(idx)
            : "";
        XYZArray.push([element, coordinates.join(" "), constraints].join(" "));
    });
    return `${XYZArray.join("\n")}\n`;
}

/**
 * Create XYZ from Material class instance (or its JSON config).
 * @param materialOrConfig {Material|Object} Material.
 * @param fractional {Boolean} Coordinate units as fractional.
 * @return {ConstrainedBasis} Class Instance
 */
function fromMaterial(materialOrConfig, fractional = false) {
    const lattice = new Lattice(materialOrConfig.lattice);
    const basis = new ConstrainedBasis({
        ...materialOrConfig.basis,
        cell: lattice.vectorArrays,
    });
    if (fractional) {
        basis.toCrystal();
    }
    return fromBasis(basis, "%11.6f");
}

export default {
    validateFormat,
    fromMaterial,
    toBasisConfig,
    fromBasis,
    CombinatorialBasis,
    validateNumberOfAtoms,
    validateAll,
};
