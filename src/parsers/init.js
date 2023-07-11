// eslint-disable-next-line max-classes-per-file
export class BaseParser {
    constructor(options) {
        this.options = options;
    }

    // eslint-disable-next-line class-methods-use-this
    parse() {
        throw new Error("parse() is implemented in children");
    }
}

export class MaterialParser extends BaseParser {
    parse(content, property_name = "material") {
        if (!(property_name === "material")) throw new Error("Implemented for material only");
        return this.parseMaterial(content);
    }

    parseMaterial(content) {
        this.content = content;
        throw new Error("parseMaterial() is implemented in children");
    }

    getCell() {
        throw new Error("Implement in children");
    }

    getElements() {
        throw new Error("Implement in children");
    }

    getCoordinates() {
        throw new Error("Implement in children");
    }

    getConstraints() {
        throw new Error("Implement in children");
    }
}
