class BaseParser {
    static property = {};

    constructor(content) {
        this.content = content;
    }

    static parse() {
        throw new Error("parse() is not defined");
    }

    /**
     * Serializes intermediate format to Object.
     * @param {Object} intermediateFormat
     * @returns {{cell: Object, elements: Object, coordinates: Object[], units: String, constraints: Object[], name: String}}
     */
    static serialize(intermediateFormat) {
        const { cell, elements, coordinates, units, constraints, name } = intermediateFormat.data;
        return { cell, elements, coordinates, units, constraints, name };
    }

    // eslint-disable-next-line class-methods-use-this
    validate() {
        throw new Error("Not Defined");
    }
}

export default BaseParser;
