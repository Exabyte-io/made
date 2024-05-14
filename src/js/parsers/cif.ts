export interface Meta {
    icsdId?: number;
}

/**
 * Extracts meta information from a CIF file
 * @param txt - CIF file text.
 */
function parseMeta(txt: string): Meta {
    const REGEX = {
        ICSD_ID: /_database_code_ICSD\s+(\d+)/,
    };
    const meta: Meta = {};
    const groups = txt.match(REGEX.ICSD_ID);

    if (groups) {
        meta.icsdId = Number(groups[1]);
    }

    return meta;
}

export default {
    parseMeta,
};
