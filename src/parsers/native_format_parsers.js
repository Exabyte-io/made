import { ConstrainedBasis } from "../basis/constrained_basis";
import { Lattice } from "../lattice/lattice";
import EspressoParser from "./espresso/init";
import Poscar from "./poscar";

const NATIVE_FORMAT = {
    JSON: "json",
    POSCAR: "poscar",
    CIF: "cif",
    QE: "qe",
    XYZ: "xyz",
    UNKNOWN: "unknown",
};

/**
 * @summary Detects the format of the input string
 * @throws {Error} - If the input string is unknown format
 * @param {String} text -  input string to detect format
 * @returns {{name: String, version:String}} - Format of the input string
 */
function detectFormat(text) {
    const jsonRegex = /^\s*\{/;
    if (jsonRegex.test(text)) return { name: NATIVE_FORMAT.JSON, version: "" };
    if (Poscar.isPoscar(text)) return { name: NATIVE_FORMAT.POSCAR, version: "5" };
    if (EspressoParser.validate(text))
        return { name: NATIVE_FORMAT.QE, version: EspressoParser.getVersionByContent(text) };
    return { name: NATIVE_FORMAT.UNKNOWN, version: "" };
}

/**
 * @summary Return intermediate format from input text
 * @param {String} text - input string to detect format and convert
 * @throws {Error} - If the input string is of unknown format
 * @return {Object} - Intermediate format
 */
function convertFromNativeFormat(text) {
    const format = detectFormat(text);

    switch (format.name) {
        case NATIVE_FORMAT.JSON:
            return JSON.parse(text);
        case NATIVE_FORMAT.POSCAR:
            return Poscar.fromPoscar(text);
        case NATIVE_FORMAT.QE:
            return EspressoParser.getIntermediateFormat(text);
        case NATIVE_FORMAT.UNKNOWN:
            throw new Error(`Unknown format`);
        // TODO:  add more formats
        default:
            throw new Error(`Unsupported format: ${format.name}`);
    }
}

function serialize(intermediateFormat) {
    // depending on metadata fetches correct parser, and its serializer to return properties
    // for now only Espresso:
    let result;
    switch (intermediateFormat.metadata.format) {
        case NATIVE_FORMAT.QE:
            result = EspressoParser.serialize(intermediateFormat);
            break;
        default:
            throw new Error(`Unsupported format: ${intermediateFormat.metadata.format}`);
    }
    const { cell, elements, coordinates, units, constraints, name } = result;
    return { cell, elements, coordinates, units, constraints, name };
}

function getMaterialConfig(text) {
    const intermediateFormat = convertFromNativeFormat(text);
    const { cell, elements, coordinates, units, constraints, name } = serialize(intermediateFormat);

    const lattice = Lattice.fromVectors({
        a: cell.vectors[0],
        b: cell.vectors[1],
        c: cell.vectors[2],
        alat: cell.alat,
        units: cell.units,
        type: cell.type,
    });

    const basis = new ConstrainedBasis({
        elements,
        coordinates,
        units,
        type: cell.type,
        cell: lattice.vectorArrays,
        constraints,
    });

    return {
        lattice: lattice.toJSON(),
        basis: basis.toJSON(),
        name,
        isNonPeriodic: false,
    };
}

export default {
    detectFormat,
    convertFromNativeFormat,
    getMaterialConfig,
};
