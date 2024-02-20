declare namespace _default {
    export { parseMeta };
}
export default _default;
/**
 * Extracts meta information from a CIF file
 * @param {String} txt - CIF file text.
 * @return {Object}
 */
declare function parseMeta(txt: string): Object;
