declare namespace _default {
    export { detectFormat };
    export { convertFromNativeFormat };
}
export default _default;
/**
 * @summary Detects the format of the input string
 * @throws {Error} - If the input string is unknown format
 * @param {string} text -  input string to detect format
 * @returns {NATIVE_FORMAT} - Format of the input string
 */
declare function detectFormat(text: string): {
    JSON: string;
    POSCAR: string;
    CIF: string;
    PWX: string;
    XYZ: string;
    UNKNOWN: string;
};
/**
 * @summary Function to handle conversion from native formats
 * @param {String} text - input string to detect format and convert
 * @throws {Error} - If the input string is of unknown format
 * @return {Object} - Material config
 */
declare function convertFromNativeFormat(text: string): Object;
