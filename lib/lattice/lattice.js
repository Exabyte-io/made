"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.Lattice = exports.nonPeriodicLatticeScalingFactor = void 0;
const lodash_1 = __importDefault(require("lodash"));
const underscore_1 = __importDefault(require("underscore"));
const cell_1 = require("../cell/cell");
const primitive_cell_1 = require("../cell/primitive_cell");
const constants_1 = require("../constants");
const math_1 = __importDefault(require("../math"));
const lattice_bravais_1 = require("./lattice_bravais");
const lattice_vectors_1 = require("./lattice_vectors");
const types_1 = require("./types");
const unit_cell_1 = require("./unit_cell");
/**
 * Scaling factor used to calculate the new lattice size for non-periodic systems.
 * The scaling factor ensures that a non-periodic structure will have have a lattice greater than the structures size.
 */
exports.nonPeriodicLatticeScalingFactor = 2.0;
/*
 * Container class for crystal lattice and associated methods.
 * Follows Bravais convention for lattice types and contains lattice vectors within.
 * Units for lattice vector coordinates are "angstroms", and "degrees" for the corresponding angles.
 */
class Lattice extends lattice_bravais_1.LatticeBravais {
    /**
     * Create a Lattice class from a config object.
     * @param {Object} config - Config object. See LatticeVectors.fromBravais.
     */
    constructor(config = {}) {
        super(config);
        this.vectors = lattice_vectors_1.LatticeVectors.fromBravais(config);
    }
    /**
     * Create a Lattice class from a list of vectors.
     * @param {Object} config - Config object. See LatticeBravais.fromVectors.
     */
    static fromVectors(config) {
        return new Lattice(lattice_bravais_1.LatticeBravais.fromVectors(config).toJSON());
    }
    /**
     * Serialize class instance to JSON.
     * @example As below:
         {
            "a" : 3.867,
            "b" : 3.867,
            "c" : 3.867,
            "alpha" : 60,
            "beta" : 60,
            "gamma" : 60,
            "units" : {
                "length" : "angstrom",
                "angle" : "degree"
            },
            "type" : "FCC",
            "vectors" : {
                "a" : [
                    3.34892,
                    0,
                    1.9335
                ],
                "b" : [
                    1.116307,
                    3.157392,
                    1.9335
                ],
                "c" : [
                    0,
                    0,
                    3.867
                ],
                "alat" : 1,
                "units" : "angstrom"
            }
        }
     */
    toJSON(skipRounding = false) {
        const round = skipRounding ? (x) => x : Lattice._roundValue; // round values by default
        return {
            a: round(this.a),
            b: round(this.b),
            c: round(this.c),
            alpha: round(this.alpha),
            beta: round(this.beta),
            gamma: round(this.gamma),
            units: {
                length: this.units.length,
                angle: this.units.angle,
            },
            type: this.type,
            vectors: this.vectors.toJSON(),
        };
    }
    clone(extraContext) {
        return new this.constructor({ ...this.toJSON(), ...extraContext });
    }
    /**
     * Get lattice vectors as a nested array
     */
    get vectorArrays() {
        return this.vectors.vectorArrays;
    }
    get Cell() {
        return new cell_1.Cell(this.vectorArrays);
    }
    /**
     * Get a short label for the type of the lattice, eg. "MCLC".
     */
    get typeLabel() {
        return lodash_1.default.get(types_1.LATTICE_TYPE_CONFIGS.find((c) => c.code === this.type), "label", "Unknown");
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
        const f_ = Lattice._roundValue;
        // construct new primitive cell using lattice parameters and skip rounding the vectors
        const pCell = (0, primitive_cell_1.primitiveCell)(latticeConfig, true);
        // create new lattice from primitive cell
        const newLattice = { ...Lattice.fromVectorArrays(pCell, latticeConfig.type) };
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
    // TODO: remove
    get unitCell() {
        const vectors = [...underscore_1.default.flatten(this.vectorArrays), this.units.length];
        return new unit_cell_1.UnitCell(vectors);
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
}
exports.Lattice = Lattice;
//# sourceMappingURL=lattice.js.map