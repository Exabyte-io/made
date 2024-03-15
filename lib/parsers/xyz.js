"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.validate = void 0;
const underscore_1 = __importDefault(require("underscore"));
const underscore_string_1 = __importDefault(require("underscore.string"));
const basis_1 = require("../basis/basis");
const constrained_basis_1 = require("../basis/constrained_basis");
const lattice_1 = require("../lattice/lattice");
const math_1 = __importDefault(require("../math"));
const errors_1 = require("./errors");
const xyz_combinatorial_basis_1 = require("./xyz_combinatorial_basis");
// Regular expression for an XYZ line with atomic constraints, eg. Si    0.000000    0.500000    0.446678 1 1 1`
// eslint-disable-next-line max-len
const XYZ_LINE_REGEX = /[A-Z][a-z]?\s+((-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)(\s+)?(\s+[0-1]\s+[0-1]\s+[0-1](\s+)?)?)$/;
/**
 * Validates XYZ file's line. Line should be in following format "Si 0.5 0.5 0.5".
 * Raises an error if line is in wrong format.
 */
function validateLine(xyzLine, index) {
    const words = xyzLine.split(" ").filter((x) => x.trim() !== "");
    if (!xyzLine.match(XYZ_LINE_REGEX)) {
        throw new errors_1.InvalidLineError(index, xyzLine);
    }
    const coordinates = [parseFloat(words[1]), parseFloat(words[2]), parseFloat(words[3])];
    coordinates.forEach((num, i) => {
        if (underscore_1.default.isNaN(num)) {
            throw new Error(`Coordinates should be a number. Possible error in ${i} coordinate`);
        }
    });
}
/**
 * Validates that passed string is well-formed XYZ file.
 */
function validate(xyzTxt) {
    // The @types/underscore.string package does not provide interfaces for function chaining at the moment. Disable TS here
    // @ts-ignore
    (0, underscore_string_1.default)(xyzTxt)
        .trim()
        .lines()
        .filter((x) => x.trim() !== "")
        .forEach(validateLine);
}
exports.validate = validate;
/**
 * Parses XYZ line and returns an object.
 * @param line - line of text
 */
function _parseXYZLineAsWords(line) {
    const words = underscore_string_1.default.words(line);
    const constraint = (n) => parseInt(`${n}`, 10) !== 0;
    return {
        element: words[0],
        coordinates: [+words[1], +words[2], +words[3]],
        // Below maps zero values to false (atom is fixed) and non-zero values to true (atom is moving)
        constraints: [constraint(+words[4]), constraint(+words[5]), constraint(+words[6])],
    };
}
/**
 * Parse XYZ text for basis. Assuming only xyz lines without blank or comment lines inbetween.
 * @param txt Text
 * @param units Coordinate units
 * @param cell Basis Cell
 */
function toBasisConfig(txt, units = "angstrom", cell = basis_1.Basis.defaultCell) {
    // @ts-ignore
    const lines = (0, underscore_string_1.default)(txt).trim().lines();
    const listOfObjects = underscore_1.default.map(lines, _parseXYZLineAsWords);
    const linePosition = (lineNumber) => lines
        .filter((l, i) => i < lineNumber)
        .map((l) => l.length + 1)
        .reduce((sum, len) => sum + len, 0);
    return {
        elements: listOfObjects.map((elm, idx) => {
            return {
                id: idx,
                value: elm.element,
                from: linePosition(idx),
                to: linePosition(idx + 1),
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
function selectionToBasis(basis, selection) {
    const aBetween = (a, b, c) => b <= a && a < c;
    basis.elements.forEach((e) => {
        e.selection = selection.ranges.some(
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion
        (r) => aBetween(r.from, e.from, e.to) || aBetween(r.to, e.from, e.to))
            ? 1
            : 0;
    });
    return basis;
}
/**
 * Create XYZ from Basis class instance.
 * @param basisClsInstance Basis class instance.
 * @param printFormat Output format for coordinates.
 * @param skipRounding Whether to round the numbers (ie. to avoid negative zeros).
 * @return Basis string in XYZ format
 */
function fromBasis(basisClsInstance, printFormat = "%9.5f", skipRounding = false) {
    const XYZArray = [];
    basisClsInstance._elements.array.forEach((item, idx) => {
        // assume that _elements and _coordinates are indexed equivalently
        const element = underscore_string_1.default.sprintf("%-3s", item);
        const coordinates = basisClsInstance.getCoordinateByIndex(idx).map((x) => {
            return underscore_string_1.default.sprintf(printFormat, skipRounding ? x : math_1.default.precise(math_1.default.roundToZero(x)));
        });
        const constraints = basisClsInstance.constraints
            ? basisClsInstance.AtomicConstraints.getAsStringByIndex(idx)
            : "";
        XYZArray.push([element, coordinates.join(" "), constraints].join(" "));
    });
    return `${XYZArray.join("\n")}\n`;
}
/**
 * Create XYZ from Material class instance (or its JSON config).
 * @param materialOrConfig Material.
 * @param fractional Coordinate units as fractional.
 * @return Class Instance
 */
function fromMaterial(materialOrConfig, fractional = false) {
    const lattice = new lattice_1.Lattice(materialOrConfig.lattice);
    // @ts-ignore
    const basis = new constrained_basis_1.ConstrainedBasis({
        ...materialOrConfig.basis,
        cell: lattice.vectorArrays,
    });
    if (fractional) {
        basis.toCrystal();
    }
    return fromBasis(basis, "%11.6f");
}
exports.default = {
    validate,
    fromMaterial,
    toBasisConfig,
    selectionToBasis,
    fromBasis,
    CombinatorialBasis: xyz_combinatorial_basis_1.CombinatorialBasis,
};
