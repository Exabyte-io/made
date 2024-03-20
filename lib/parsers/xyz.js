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
const XYZ_LINE_REGEX = /[A-Z][a-z]?(\d)?\s+((-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)(\s+)?(\s+[0-1]\s+[0-1]\s+[0-1](\s+)?)?)$/;
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
    const elementWithOptionalLabel = words[0];
    const element = elementWithOptionalLabel.replace(/\d$/, ""); // Fe1 => Fe
    const constraint = (n) => parseInt(`${n}`, 10) !== 0;
    const basisLineConfig = {
        element,
        coordinates: [+words[1], +words[2], +words[3]],
        // Below maps zero values to false (atom is fixed) and non-zero values to true (atom is moving)
        constraints: [constraint(+words[4]), constraint(+words[5]), constraint(+words[6])],
    };
    if (elementWithOptionalLabel !== element) {
        return {
            ...basisLineConfig,
            label: parseInt(elementWithOptionalLabel[elementWithOptionalLabel.length - 1], 10),
        };
    }
    return basisLineConfig;
}
/**
 * Parse XYZ text for basis.
 * @param txt Text
 * @param units Coordinate units
 * @param cell Basis Cell
 */
function toBasisConfig(txt, units = "angstrom", cell = basis_1.Basis.defaultCell) {
    // @ts-ignore
    const lines = (0, underscore_string_1.default)(txt).trim().lines();
    const listOfObjects = underscore_1.default.map(lines, _parseXYZLineAsWords);
    const basisConfig = {
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
    const labels = [];
    listOfObjects.forEach((elm, idx) => {
        if (elm.label) {
            labels.push({
                id: idx,
                value: elm.label,
            });
        }
    });
    if (!underscore_1.default.isEmpty(labels)) {
        return {
            ...basisConfig,
            labels,
        };
    }
    return basisConfig;
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
        const atomicLabel = basisClsInstance.atomicLabelsArray[idx];
        const elementWithLabel = item + atomicLabel;
        const element = underscore_string_1.default.sprintf("%-3s", elementWithLabel);
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
    fromBasis,
    CombinatorialBasis: xyz_combinatorial_basis_1.CombinatorialBasis,
};
