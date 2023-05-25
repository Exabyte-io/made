/**
 * @summary Function to handle conversion from native formats
 * @param {String} text - Text to autodetect format and convert
 * @throws {Error} - If the input is not a valid JSON string
 * @return {Object} - Material config
 */
function convertFromNative(text) {
    let materialConfig = {};
    try {
        materialConfig = JSON.parse(text);
    } catch (error) {
        throw new Error("Invalid JSON: " + error.message);
    }
    return materialConfig;
}

export default convertFromNative;
