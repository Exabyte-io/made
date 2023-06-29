import { APPLICATIONS, NATIVE_FORMATS } from "../enums";
import BaseParser from "../init";
import parser7_2 from "./7.2/parser";

const parserMap = {
    // TODO: change correctly
    5.2: parser7_2,
    "5.*": parser7_2,
    "6.*": parser7_2,
    7.2: parser7_2,
};

// eslint-disable-next-line no-unused-vars
const versionRegexMap = {
    // TODO: implement version detection
    5.2: /temp/,
    "5.*": function () {},
    "6.*": function () {},
    "7.*": function () {},
};

class EspressoParser extends BaseParser {
    // eslint-disable-next-line no-useless-constructor
    constructor(content) {
        super(content);
    }

    static validate(text) {
        const version = EspressoParser.getVersionByContent(text);
        const parser = EspressoParser.getParserByVersion(version);
        return parser.isEspressoFormat(text);
    }

    // eslint-disable-next-line no-unused-vars
    static getVersionByContent(text) {
        // eslint-disable-next-line no-restricted-syntax
        // for (const version in versionRegexMap) {
        //     if (versionRegexMap[version].test(text)) {
        //         return version;
        //     }
        // }
        // throw new Error("Could not determine version from content.");
        return "7.2"; // TODO: REMOVE
    }

    static getParserByVersion(version) {
        return parserMap[version];
    }

    static getIntermediateFormat(text) {
        const version = EspressoParser.getVersionByContent(text);
        const parser = EspressoParser.getParserByVersion(version);
        return {
            metadata: {
                application: APPLICATIONS.ESPRESSO,
                format: NATIVE_FORMATS.QE,
                version: EspressoParser.getVersionByContent(text),
            },
            data: parser.fromEspressoFormat(text),
        };
    }

    static serialize(intermediateFormat) {
        // get all the paramteres according to the version
        const { cell, elements, coordinates, units, constraints, name } = intermediateFormat.data;
        return { cell, elements, coordinates, units, constraints, name };
    }
}

export default EspressoParser;
