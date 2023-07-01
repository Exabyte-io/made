import { ConstrainedBasis } from "../basis/constrained_basis";
import { Lattice } from "../lattice/lattice";
import { NATIVE_FORMATS } from "./enums";
import EspressoParser from "./espresso/init";
import Poscar from "./poscar";
import PoscarParser from "./poscar/init";

/**
 * @summary Detects the format of the input string
 * @throws {Error} - If the input string is unknown format
 * @param {String} text -  input string to detect format
 * @returns {{name: String, version:String}} - Format of the input string
 */
function detectFormat(text) {
    const jsonRegex = /^\s*\{/;
    if (jsonRegex.test(text)) return { name: NATIVE_FORMATS.JSON, version: "" };
    if (PoscarParser.validate(text)) return { name: NATIVE_FORMATS.POSCAR, version: "5" };
    if (EspressoParser.validate(text))
        return { name: NATIVE_FORMATS.QE, version: EspressoParser.getVersionByContent(text) };
    return { name: NATIVE_FORMATS.UNKNOWN, version: "" };
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
        case NATIVE_FORMATS.JSON:
            return JSON.parse(text);
        case NATIVE_FORMATS.POSCAR:
            return Poscar.fromPoscar(text); // TODO: rewrite when ready to use PoscarParser
        // return PoscarParser.getIntermediateFormat(text);
        case NATIVE_FORMATS.QE:
            return EspressoParser.getIntermediateFormat(text);
        case NATIVE_FORMATS.UNKNOWN:
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
        case NATIVE_FORMATS.QE:
            result = EspressoParser.serialize(intermediateFormat);
            break;
        case NATIVE_FORMATS.POSCAR:
            result = PoscarParser.serialize(intermediateFormat);
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
