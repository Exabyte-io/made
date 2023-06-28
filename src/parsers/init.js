class BaseParser {
    static property = {};

    constructor(content) {
        this.content = content;
    }

    // Returns desired object
    // eslint-disable-next-line class-methods-use-this
    serialize() {
        throw new Error("Not Defined");
    }

    // eslint-disable-next-line class-methods-use-this
    validate() {
        throw new Error("Not Defined");
    }
}

export default BaseParser;
