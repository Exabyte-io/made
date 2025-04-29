"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.ElementWithLabel = void 0;
const underscore_string_1 = __importDefault(require("underscore.string"));
class ElementWithLabel {
    constructor({ element, label }) {
        this.element = element;
        this.label = label;
    }
    prettyPrint() {
        const elementWithLabel = `${this.element}${this.label}`;
        return underscore_string_1.default.sprintf("%-3s", elementWithLabel);
    }
}
exports.ElementWithLabel = ElementWithLabel;
