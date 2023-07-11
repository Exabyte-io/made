import { STRUCTURAL_INFORMATION_FORMATS } from "./enums";
import { ESPRESSOMaterialParser } from "./espresso/parser";
import Poscar from "./poscar";
/**
 * @summary Detects the format of the input string
 * @throws {Error} - If the input string is unknown format
 * @param {string} text -  input string to detect format
 * @returns {string} - Format of the input string
 */
function detectFormat(text) {
    const jsonRegex = /^\s*\{/;
    const espressoRegex = /^\s*ATOMIC_SPECIES/; // TODO: replace with actual detection function
    if (jsonRegex.test(text)) return STRUCTURAL_INFORMATION_FORMATS.JSON;
    if (Poscar.isPoscar(text)) return STRUCTURAL_INFORMATION_FORMATS.POSCAR;
    if (espressoRegex.test(text)) return STRUCTURAL_INFORMATION_FORMATS.QE;
    return STRUCTURAL_INFORMATION_FORMATS.UNKNOWN;
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
        case STRUCTURAL_INFORMATION_FORMATS.JSON:
            return JSON.parse(text);
        case STRUCTURAL_INFORMATION_FORMATS.POSCAR:
            return Poscar.fromPoscar(text);
        case STRUCTURAL_INFORMATION_FORMATS.QE:
            // eslint-disable-next-line no-case-declarations
            const parser = new ESPRESSOMaterialParser(); // TODO: replace with parsers factory
            return parser.parse(text, "material");
        case STRUCTURAL_INFORMATION_FORMATS.UNKNOWN:
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
