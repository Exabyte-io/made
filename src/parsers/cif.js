/**
 * @summary Extracts cif meta information
 */
function parseMeta(txt) {
    const REGEX = {
        ICSD_ID: /_database_code_ICSD\s+(\d+)/
    };

    const meta = {};
    const groups = txt.match(REGEX.ICSD_ID);

    if (groups) {
        meta.icsdId = Number(groups[1]);
    }

    return meta;
}

export default {
    parseMeta
};
