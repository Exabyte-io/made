/**
 * @summary Replaces placeholder instances within a string with specified values.
 * The placeholders are denoted as {n}, where n is the index of the value to be replaced in the array or arguments passed.
 * If a placeholder has no corresponding value, the placeholder remains in the string.
 * @param {String} str
 * @param {String | String[]} values
 * @returns {String}
 */
export function formatString(str, ...values) {
    return str.replace(/{(\d+)}/g, (match, index) => {
        return typeof values[index] !== "undefined" ? values[index] : match;
    });
}
