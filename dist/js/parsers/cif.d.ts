export interface Meta {
    icsdId?: number;
}
/**
 * Extracts meta information from a CIF file
 * @param txt - CIF file text.
 */
declare function parseMeta(txt: string): Meta;
declare const _default: {
    parseMeta: typeof parseMeta;
};
export default _default;
