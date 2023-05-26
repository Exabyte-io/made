import Poscar from "./poscar";

// TODO: move this enum to an appropriate file
export const NATIVE_FORMAT = {
    JSON: "json",
    POSCAR: "poscar",
    CIF: "cif",
    PWX: "pwx",
    XYZ: "xyz",
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

    throw new Error("Unknown format");
}

/**
 * @summary Function to handle conversion from native formats
 * @param {String} text - input string to detect format and convert
 * @throws {Error} - If the input string is of unknown format
 * @return {Object} - Material config
 */
function convertFromNative(text) {
    const format = detectFormat(text);

    switch (format) {
        case NATIVE_FORMAT.JSON:
            return JSON.parse(text);
        case NATIVE_FORMAT.POSCAR:
            return Poscar.fromPoscar(text);
        // TODO:  add more formats
        default:
            throw new Error(`Unsupported format: ${format}`);
    }
}

export default convertFromNative;
