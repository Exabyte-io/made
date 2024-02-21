export class BaseParser {
    constructor(options) {
        this.options = options;
    }

    // eslint-disable-next-line class-methods-use-this
    parse() {
        throw new Error("parse() is implemented in children");
    }
}
