"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
const cif_1 = __importDefault(require("./cif"));
const espresso_1 = __importDefault(require("./espresso"));
const native_format_parsers_1 = __importDefault(require("./native_format_parsers"));
const poscar_1 = __importDefault(require("./poscar"));
const xyz_1 = __importDefault(require("./xyz"));
exports.default = {
    xyz: xyz_1.default,
    poscar: poscar_1.default,
    cif: cif_1.default,
    espresso: espresso_1.default,
    nativeFormatParsers: native_format_parsers_1.default,
};
//# sourceMappingURL=parsers.js.map