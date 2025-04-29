"use strict";
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    var desc = Object.getOwnPropertyDescriptor(m, k);
    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
      desc = { enumerable: true, get: function() { return m[k]; } };
    }
    Object.defineProperty(o, k2, desc);
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar = (this && this.__importStar) || function (mod) {
    if (mod && mod.__esModule) return mod;
    var result = {};
    if (mod != null) for (var k in mod) if (k !== "default" && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
    __setModuleDefault(result, mod);
    return result;
};
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.Lattice = exports.LatticeVectors = exports.nonPeriodicLatticeScalingFactor = void 0;
const constants_1 = require("@mat3ra/code/dist/js/constants");
const entity_1 = require("@mat3ra/code/dist/js/entity");
const lodash = __importStar(require("lodash"));
const cell_1 = require("../cell/cell");
const primitive_cell_1 = require("../cell/primitive_cell");
const math_1 = __importDefault(require("../math"));
const lattice_types_1 = require("./lattice_types");
const unit_cell_1 = require("./unit_cell");
/**
 * Scaling factor used to calculate the new lattice size for non-periodic systems.
 * The scaling factor ensures that a non-periodic structure will have have a lattice greater than the structures size.
 */
exports.nonPeriodicLatticeScalingFactor = 2.0;
class LatticeVectors extends cell_1.Cell {
}
exports.LatticeVectors = LatticeVectors;
class Lattice extends entity_1.InMemoryEntity {
    constructor(config = Lattice.defaultConfig) {
        super(config);
        this.type = "TRI";
        const { a, b, c, alpha, beta, gamma, units, type } = config;
        this.a = a;
        this.b = b;
        this.c = c;
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;
        this.units = units;
        this.type = type || "TRI";
    }
    static fromConfig(config) {
        const latticeConfig = config;
        return new Lattice(latticeConfig);
    }
    static fromConfigPartial(config) {
        const primitiveLatticeConfig = Lattice.getDefaultPrimitiveLatticeConfigByType(config);
        return new Lattice(primitiveLatticeConfig);
    }
    calculateVectors() {
        const { a } = this;
        const { b } = this;
        const { c } = this;
        const alphaRad = math_1.default.unit(this.alpha, "deg").toNumber("rad");
        const betaRad = math_1.default.unit(this.beta, "deg").toNumber("rad");
        const gammaRad = math_1.default.unit(this.gamma, "deg").toNumber("rad");
        const cosAlpha = math_1.default.cos(alphaRad);
        const cosBeta = math_1.default.cos(betaRad);
        const cosGamma = math_1.default.cos(gammaRad);
        const sinAlpha = math_1.default.sin(alphaRad);
        const sinBeta = math_1.default.sin(betaRad);
        const gammaStar = math_1.default.acos((cosAlpha * cosBeta - cosGamma) / (sinAlpha * sinBeta));
        const cosGammaStar = math_1.default.cos(gammaStar);
        const sinGammaStar = math_1.default.sin(gammaStar);
        const vectorA = [a * sinBeta, 0.0, a * cosBeta];
        const vectorB = [
            -b * sinAlpha * cosGammaStar,
            b * sinAlpha * sinGammaStar,
            b * cosAlpha,
        ];
        const vectorC = [0.0, 0.0, c];
        return [vectorA, vectorB, vectorC];
    }
    static fromVectors(config) {
        return this.fromVectorsArray([config.a, config.b, config.c]);
    }
    static fromVectorsArray(vectors, units = Lattice.defaultConfig.units, type = "TRI") {
        const [aVec, bVec, cVec] = vectors;
        const a = math_1.default.vlen(aVec);
        const b = math_1.default.vlen(bVec);
        const c = math_1.default.vlen(cVec);
        const alpha = math_1.default.angle(bVec, cVec, "deg");
        const beta = math_1.default.angle(aVec, cVec, "deg");
        const gamma = math_1.default.angle(aVec, bVec, "deg");
        return new Lattice({
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
            units,
            type,
        });
    }
    // TODO: remove
    get unitCell() {
        var _a;
        const units = ((_a = this.units) === null || _a === void 0 ? void 0 : _a.length) || "angstrom";
        const vectors = [...lodash.flatten(this.vectorArrays), units];
        return new unit_cell_1.UnitCell(vectors);
    }
    get vectors() {
        return LatticeVectors.fromVectorsArray(this.calculateVectors());
    }
    get vectorArrays() {
        return this.vectors.vectorArrays;
    }
    get vectorArraysRounded() {
        return this.vectors.vectorArraysRounded;
    }
    get cellVolume() {
        return this.vectors.volume;
    }
    get cellVolumeRounded() {
        return this.vectors.volumeRounded;
    }
    /**
     * Get a short label for the type of the lattice, eg. "MCLC".
     */
    get typeLabel() {
        return lodash.get(lattice_types_1.LATTICE_TYPE_CONFIGS.find((c) => c.code === this.type), "label", "Unknown");
    }
    /**
     * Get a short label for the extended type of the lattice, eg. "MCLC-5".
     */
    get typeExtended() {
        const { a, b, c, alpha, beta, gamma, type } = this;
        const cosAlpha = math_1.default.cos((alpha / 180) * math_1.default.PI);
        switch (type) {
            case "BCT":
                return c < a ? "BCT-1" : "BCT-2";
            case "ORCF":
                if (1 / (a * a) >= 1 / (b * b) + 1 / (c * c)) {
                    return "ORCF-1";
                }
                return "ORCF-2";
            case "RHL":
                return cosAlpha > 0 ? "RHL-1" : "RHL-2";
            case "MCLC":
                if (gamma >= 90) {
                    // MCLC-1,2
                    return "MCLC-1";
                }
                if ((b / c) * cosAlpha + ((b * b) / (a * a)) * (1 - cosAlpha * cosAlpha) <= 1) {
                    // MCLC-3,4
                    return "MCLC-3";
                }
                return "MCLC-5";
            case "TRI":
                if (alpha > 90 && beta > 90 && gamma >= 90) {
                    // TRI-1a,2a
                    return "TRI_1a";
                }
                return "TRI_1b";
            default:
                return type;
        }
    }
    /**
     * Calculate the volume of the lattice cell.
     */
    get volume() {
        return math_1.default.abs(math_1.default.det(this.vectorArrays));
    }
    /*
     * Returns a "default" primitive lattice by type, with lattice parameters scaled by the length of "a",
     * @param latticeConfig {Object} LatticeBravais config (see constructor)
     */
    static getDefaultPrimitiveLatticeConfigByType(latticeConfig) {
        const f_ = math_1.default.roundArrayOrNumber;
        // construct new primitive cell using lattice parameters and skip rounding the vectors
        const vectors = (0, primitive_cell_1.getPrimitiveLatticeVectorsFromConfig)(latticeConfig);
        // create new lattice from primitive cell
        const newLattice = Lattice.fromVectorsArray(vectors, undefined, latticeConfig.type);
        // preserve the new type and scale back the lattice parameters
        const k = latticeConfig.a / newLattice.a;
        return Object.assign(newLattice, {
            a: f_(newLattice.a * k),
            b: f_(newLattice.b * k),
            c: f_(newLattice.c * k),
            alpha: f_(newLattice.alpha),
            beta: f_(newLattice.beta),
            gamma: f_(newLattice.gamma),
        });
    }
    /**
     * Returns a string further used for the calculation of an unique hash.
     * @param isScaled - Whether to scale the vectors by the length of the first vector initially.
     */
    getHashString(isScaled = false) {
        // lattice vectors must be measured in angstroms
        const scaleK = isScaled ? this.a : 1;
        const scaledLattice = {
            ...this,
            a: this.a / scaleK,
            b: this.b / scaleK,
            c: this.c / scaleK,
        };
        // form lattice string
        return `${[
            scaledLattice.a,
            scaledLattice.b,
            scaledLattice.c,
            scaledLattice.alpha,
            scaledLattice.beta,
            scaledLattice.gamma,
        ]
            .map((x) => math_1.default.round(x, constants_1.HASH_TOLERANCE))
            .join(";")};`;
    }
    /**
     * Get the list of editable keys (eg. 'a', 'alpha') for the current lattice.
     * @return {Object}
     * @example {a: true, b: false, c: false, alpha: true, beta: false, gamma: false}
     */
    get editables() {
        var _a;
        const object = {};
        const editablesList = (_a = lattice_types_1.LATTICE_TYPE_CONFIGS.find((entry) => entry.code === this.type)) === null || _a === void 0 ? void 0 : _a.editables;
        // ["a", "gamma"] => {a: true, gamma: true}
        if (editablesList) {
            editablesList.forEach((element) => {
                Object.assign(object, {
                    [element]: true,
                });
            });
        }
        return object;
    }
    // @ts-ignore
    toJSON(exclude) {
        return {
            ...super.toJSON(exclude),
            vectors: this.vectors.toJSON(),
        };
    }
}
exports.Lattice = Lattice;
Lattice.defaultConfig = {
    a: 1,
    b: 1,
    c: 1,
    alpha: 90,
    beta: 90,
    gamma: 90,
    units: { length: "angstrom", angle: "degree" },
    type: "CUB",
};
