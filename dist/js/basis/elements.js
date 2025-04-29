"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Elements = exports.Element = void 0;
const code_1 = require("@mat3ra/code");
const lodash_1 = require("lodash");
class Element extends code_1.ValueWithId {
    constructor({ value, id }) {
        super({ id, value });
        this.value = value;
    }
    prettyPrint() {
        return this.value;
    }
}
exports.Element = Element;
class Elements extends code_1.ArrayWithIds {
    getUnique() {
        return (0, lodash_1.uniq)(this.values.map((element) => element));
    }
}
exports.Elements = Elements;
