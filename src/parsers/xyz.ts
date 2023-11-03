import _ from "underscore";
import s from "underscore.string";

import { Basis } from "../basis/basis";
import { ConstrainedBasis } from "../basis/constrained_basis";
import { Lattice } from "../lattice/lattice";
import math from "../math";
import { InvalidLineError } from "./errors";
import { CombinatorialBasis } from "./xyz_combinatorial_basis";

// Regular expression for an XYZ line with atomic constraints, eg. Si    0.000000    0.500000    0.446678 1 1 1`
// eslint-disable-next-line max-len
const XYZ_LINE_REGEX =
    /[A-Z][a-z]?\s+((-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)(\s+)?(\s+[0-1]\s+[0-1]\s+[0-1](\s+)?)?)$/;

/**
 * Validates XYZ file's line. Line should be in following format "Si 0.5 0.5 0.5".
 * Raises an error if line is in wrong format.
 */
function validateLine(xyzLine: string, index: number) {
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
 * Validates that passed string is well-formed XYZ file.
 */
export function validate(xyzTxt: string) {
    // The @types/underscore.string package does not provide interfaces for function chaining at the moment. Disable TS here
    // @ts-ignore
    s(xyzTxt)
        .trim()
        .lines()
        .filter((x: string) => x.trim() !== "")
        .forEach(validateLine);
}

interface ParsedObject {
    element: string;
    coordinates: [number, number, number];
    constraints: [boolean, boolean, boolean];
}

/**
 * Parses XYZ line and returns an object.
 * @param line - line of text
 */
function _parseXYZLineAsWords(line: string): ParsedObject {
    const words = s.words(line);
    const constraint = (n: number) => parseInt(`${n}`, 10) !== 0;
    return {
        element: words[0],
        coordinates: [+words[1], +words[2], +words[3]],
        // Below maps zero values to false (atom is fixed) and non-zero values to true (atom is moving)
        constraints: [constraint(+words[4]), constraint(+words[5]), constraint(+words[6])],
    };
}

/**
 * Parse XYZ text for basis.
 * @param txt Text
 * @param units Coordinate units
 * @param cell Basis Cell
 */
function toBasisConfig(txt: string, units = "angstrom", cell = Basis.defaultCell) {
    // @ts-ignore
    const lines: string[] = s(txt).trim().lines();
    const listOfObjects = _.map(lines, _parseXYZLineAsWords);

    return {
        elements: listOfObjects.map((elm, idx) => {
            return {
                id: idx,
                value: elm.element,
            };
        }),
        coordinates: listOfObjects.map((elm, idx) => {
            return {
                id: idx,
                value: elm.coordinates,
            };
        }),
        units,
        cell,
        constraints: listOfObjects.map((elm, idx) => {
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
        const coordinates = clsInstance.getCoordinateByIndex(idx).map((x) => {
            return s.sprintf(printFormat, skipRounding ? x : math.precise(math.roundToZero(x)));
        });
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
 * @return Class Instance
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
    validate,
    fromMaterial,
    toBasisConfig,
    fromBasis,
    CombinatorialBasis,
};
