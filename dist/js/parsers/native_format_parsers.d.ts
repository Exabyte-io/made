/**
 * @summary Detects the format of the input string
 * @throws {Error} - If the input string is unknown format
 * @param {string} text -  input string to detect format
 * @returns {NATIVE_FORMAT} - Format of the input string
 */
declare function detectFormat(text: string): string;
/**
 * @summary Function to handle conversion from native formats
 * @param {String} text - input string to detect format and convert
 * @throws {Error} - If the input string is of unknown format
 * @return {Object} - Material config
 */
declare function convertFromNativeFormat(text: string): any;
declare const _default: {
    detectFormat: typeof detectFormat;
    convertFromNativeFormat: typeof convertFromNativeFormat;
};
export default _default;
