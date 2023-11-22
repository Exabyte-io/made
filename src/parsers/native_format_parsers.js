import Poscar from "./poscar";

const NATIVE_FORMAT = {
    JSON: "json",
    POSCAR: "poscar",
    CIF: "cif",
    PWX: "pwx",
    XYZ: "xyz",
    UNKNOWN: "unknown",
};

/**
 * @summary Detects the format of the input string
 * @throws {Error} - If the input string is unknown format
 * @param {string} text -  input string to detect format
 * @returns {NATIVE_FORMAT} - Format of the input string
 */
function detectFormat(text) {
    const jsonRegex = /^\s*\{/;
    if (jsonRegex.test(text)) return NATIVE_FORMAT.JSON;
    if (Poscar.isPoscar(text)) return NATIVE_FORMAT.POSCAR;

    return NATIVE_FORMAT.UNKNOWN;
}

/**
 * @summary Function to handle conversion from native formats
 * @param {String} text - input string to detect format and convert
 * @throws {Error} - If the input string is of unknown format
 * @return {Object} - Material config
 */
function convertFromNativeFormat(text) {
    const format = detectFormat(text);

    switch (format) {
        case NATIVE_FORMAT.JSON:
            return JSON.parse(text);
        case NATIVE_FORMAT.POSCAR:
            return Poscar.fromPoscar(text);
        case NATIVE_FORMAT.UNKNOWN:
            throw new Error(`Unknown format`);
        // TODO:  add more formats
        default:
            throw new Error(`Unsupported format: ${format}`);
    }
}

export default {
    detectFormat,
    convertFromNativeFormat,
};
