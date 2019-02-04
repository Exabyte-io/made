import _ from "underscore";
import s from "underscore.string";

import {Basis} from "../basis/basis";
import {InvalidLineError} from "./errors";
import {Lattice} from "../lattice/lattice";
import {ConstrainedBasis} from "../basis/constrained_basis";
import {CombinatorialBasis} from "./xyz_combinatorial_basis";

// Regular expression for an XYZ line with atomic constraints, eg. Si    0.000000    0.500000    0.446678 1 1 1`
const XYZ_LINE_REGEX = /[A-Z][a-z]?\s+((-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)\s+(-?\d+\.?\d*|\.\d+)(\s+[0-1]\s+[0-1]\s+[0-1]\s+?)?)/;

/**
 * @summary Validates that passed string is well-formed XYZ file.
 * @param xyzTxt {String}
 */
export function validate(xyzTxt) {
    s(xyzTxt).trim().lines().filter(x => x.trim() !== '').forEach(validateLine);
}

/**
 * @summary Validates XYZ file's line. Line should be in following format "Si 0.5 0.5 0.5".
 * Raises an error if line is in wrong format.
 * @param xyzLine {String}
 * @param index {Number}
 */
function validateLine(xyzLine, index) {
    const words = xyzLine.split(' ').filter(x => x.trim() !== '');

    if (!xyzLine.match(XYZ_LINE_REGEX)) {
        throw new InvalidLineError(index, xyzLine);
    }

    const coordinates = [parseFloat(words[1]), parseFloat(words[2]), parseFloat(words[3])];
    coordinates.forEach(function (num, index) {
        if (_.isNaN(num)) {
            throw new Error(`Coordinates should be a number. Possible error in ${i} coordinate`);
        }
    });
}

function _parseXYZLineAsWords(line) {
    const words = s.words(line);
    return {
        element: words[0],
        coordinates: [(+words[1]), (+words[2]), (+words[3])],
        constraints: [(+words[4]), (+words[5]), (+words[6])].map(e => !!e),
    };
}

/**
 * @summary Parse XYZ text for basis.
 * @param txt {String} Text
 * @param units {String} Coordinate units
 * @param cell {Array} Basis Cell
 * @return {Object}
 */
function toBasisConfig(txt, units = 'angstrom', cell = Basis.defaultCell) {
    const lines = s(txt).trim().lines();
    const listOfObjects = _.map(lines, _parseXYZLineAsWords);

    const addSelectiveDynamics = listOfObjects.some(line => (line.constraints || []).some(e => e));

    return {
        // using concat below to avoid modifying listOfObjects in place with map
        elements: [].concat(listOfObjects).map((elm, idx) => {
            return {
                id: idx,
                value: elm.element
            }
        }),
        coordinates: [].concat(listOfObjects).map((elm, idx) => {
            return {
                id: idx,
                value: elm.coordinates
            }
        }),
        units: units,
        cell: cell,
        constraints: !addSelectiveDynamics ? [] : [].concat(listOfObjects).map((elm, idx) => {
            return {
                id: idx,
                value: elm.constraints
            }
        })

    };
}

/**
 * @summary create XYZ from Basis class instance.
 * @param basisClsInstance {ConstrainedBasis} Basis class instance.
 * @param printFormat {String} Output format for coordinates.
 * @return {String} Basis string in XYZ format
 */
function fromBasis(basisClsInstance, printFormat = '%9.5f') {
    const clsInstance = basisClsInstance;
    const XYZArray = [];
    clsInstance._elements.array.forEach((item, idx) => {
        // assume that _elements and _coordinates are indexed equivalently
        const element = s.sprintf('%-3s', item);
        const coordinates = clsInstance.getCoordinateByIndex(idx).map(x => s.sprintf(printFormat, x));
        const constraints = clsInstance.constraints ? clsInstance.AtomicConstraints.getAsStringByIndex(idx) : '';
        XYZArray.push([element, coordinates.join(' '), constraints].join(' '));
    });
    return XYZArray.join('\n') + '\n';
}

/**
 * @summary create XYZ from Material class instance (or its JSON config).
 * @param materialOrConfig {Material|Object} Material.
 * @param fractional {Boolean} Coordinate units as fractional.
 * @return {ConstrainedBasis} Class Instance
 */
function fromMaterial(materialOrConfig, fractional = false) {
    const lattice = new Lattice(materialOrConfig.lattice);
    const basis = new ConstrainedBasis({
        ...materialOrConfig.basis,
        cell: lattice.vectorArrays
    });
    fractional && basis.toCrystal();
    return fromBasis(basis, '%11.6f')
}

export default {
    validate,
    fromMaterial,
    toBasisConfig,
    fromBasis,
    CombinatorialBasis,
};
