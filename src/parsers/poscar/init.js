import { APPLICATIONS, NATIVE_FORMATS } from "../enums";
import BaseParser from "../init";
import parser5 from "./5/parser";

const parserMap = {
    5: parser5,
};

// eslint-disable-next-line no-unused-vars
const versionRegexMap = {
    5: /temp/,
};

class PoscarParser extends BaseParser {
    // eslint-disable-next-line no-useless-constructor
    constructor(content) {
        super(content);
    }

    static validate(text) {
        const version = PoscarParser.getVersionByContent(text);
        const parser = PoscarParser.getParserByVersion(version);
        return parser.isPoscar(text);
    }

    // eslint-disable-next-line no-unused-vars
    static getVersionByContent(text) {
        return "5";
    }

    static getParserByVersion(version) {
        return parserMap[version];
    }

    static getIntermediateFormat(text) {
        const version = PoscarParser.getVersionByContent(text);
        const parser = PoscarParser.getParserByVersion(version);
        return {
            metadata: {
                application: APPLICATIONS.VASP,
                format: NATIVE_FORMATS.POSCAR,
                version: PoscarParser.getVersionByContent(text),
            },
            data: parser.fromPoscar(text),
        };
    }
}

export default PoscarParser;
