"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.Material = exports.MaterialMixin = exports.defaultMaterialConfig = void 0;
const entity_1 = require("@mat3ra/code/dist/js/entity");
const crypto_js_1 = __importDefault(require("crypto-js"));
const constrained_basis_1 = require("./basis/constrained_basis");
const conventional_cell_1 = require("./cell/conventional_cell");
const constants_1 = require("./constants");
const lattice_1 = require("./lattice/lattice");
const parsers_1 = __importDefault(require("./parsers/parsers"));
// TODO: fix dependency cycle below
// eslint-disable-next-line import/no-cycle
const supercell_1 = __importDefault(require("./tools/supercell"));
exports.defaultMaterialConfig = {
    name: "Silicon FCC",
    basis: {
        elements: [
            {
                id: 1,
                value: "Si",
            },
            {
                id: 2,
                value: "Si",
            },
        ],
        coordinates: [
            {
                id: 1,
                value: [0.0, 0.0, 0.0],
            },
            {
                id: 2,
                value: [0.25, 0.25, 0.25],
            },
        ],
        units: constants_1.ATOMIC_COORD_UNITS.crystal,
    },
    lattice: {
        // Primitive cell for Diamond FCC Silicon at ambient conditions
        type: "FCC",
        a: 3.867,
        b: 3.867,
        c: 3.867,
        alpha: 60,
        beta: 60,
        gamma: 60,
        units: {
            length: constants_1.units.angstrom,
            angle: constants_1.units.degree,
        },
    },
};
function MaterialMixin(superclass) {
    class MadeMaterial extends superclass {
        // TODO: add constraints (and other properties if needed) to ESSE MaterialSchema, then uncomment the line below to allow validation
        // During validation of the Material entity, properties absent in ESSE schema get deleted.
        // static readonly jsonSchema = MaterialJSONSchemaObject;
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        constructor(...config) {
            super(...config);
            this.name = super.name || this.formula;
        }
        toJSON() {
            return {
                ...super.toJSON(),
                lattice: this.Lattice.toJSON(),
                basis: this.Basis.toJSON(),
                name: this.name,
                isNonPeriodic: this.isNonPeriodic,
            };
        }
        static get defaultConfig() {
            return exports.defaultMaterialConfig;
        }
        get src() {
            return this.prop("src");
        }
        set src(src) {
            this.setProp("src", src);
        }
        updateFormula() {
            this.setProp("formula", this.Basis.formula);
            this.setProp("unitCellFormula", this.Basis.unitCellFormula);
        }
        /**
         * Gets Bolean value for whether or not a material is non-periodic vs periodic.
         * False = periodic, True = non-periodic
         */
        get isNonPeriodic() {
            return this.prop("isNonPeriodic", false);
        }
        /**
         * @summary Sets the value of isNonPeriodic based on Boolean value passed as an argument.
         */
        set isNonPeriodic(bool) {
            this.setProp("isNonPeriodic", bool);
        }
        /**
         * @summary Returns the specific derived property (as specified by name) for a material.
         */
        getDerivedPropertyByName(name) {
            return this.getDerivedProperties().find((x) => x.name === name);
        }
        /**
         * @summary Returns the derived properties array for a material.
         */
        getDerivedProperties() {
            return this.prop("derivedProperties", []);
        }
        /**
         * Gets material's formula
         */
        get formula() {
            return this.prop("formula") || this.Basis.formula;
        }
        get unitCellFormula() {
            return this.prop("unitCellFormula") || this.Basis.unitCellFormula;
        }
        // should be private, but TS throws error "Property 'unsetFileProps' of exported class expression may not be private or protected"
        unsetFileProps() {
            this.unsetProp("src");
            this.unsetProp("icsdId");
            this.unsetProp("external");
        }
        /**
         * @param textOrObject Basis text or JSON object.
         * @param format Format (xyz, etc.)
         * @param unitz crystal/cartesian
         */
        setBasis(textOrObject, format, unitz) {
            let basis;
            switch (format) {
                case "xyz":
                    basis = parsers_1.default.xyz.toBasisConfig(textOrObject, unitz);
                    break;
                default:
                    basis = textOrObject;
            }
            this.setProp("basis", basis);
            this.unsetFileProps();
            this.updateFormula();
        }
        setBasisConstraints(constraints) {
            this.setBasis({ ...this.basis, constraints });
        }
        get basis() {
            return this.prop("basis");
        }
        // returns the instance of {ConstrainedBasis} class
        get Basis() {
            return new constrained_basis_1.ConstrainedBasis({
                ...this.basis,
                cell: this.Lattice.vectorArrays,
            });
        }
        /**
         * High-level access to unique elements from material instead of basis.
         */
        get uniqueElements() {
            return this.Basis.uniqueElements;
        }
        get lattice() {
            return this.prop("lattice", undefined);
        }
        set lattice(config) {
            this.setProp("lattice", config);
            this.unsetFileProps();
        }
        get Lattice() {
            return new lattice_1.Lattice(this.lattice);
        }
        /**
         * Returns the inchi string from the derivedProperties for a non-periodic material, or throws an error if the
         *  inchi cannot be found.
         *  @returns {String}
         */
        getInchiStringForHash() {
            const inchi = this.getDerivedPropertyByName("inchi");
            if (inchi) {
                return inchi.value;
            }
            throw new Error("Hash cannot be created. Missing InChI string in derivedProperties");
        }
        /**
         * Calculates hash from basis and lattice. Algorithm expects the following:
         * - asserts lattice units to be angstrom
         * - asserts basis units to be crystal
         * - asserts basis coordinates and lattice measurements are rounded to hash precision
         * - forms strings for lattice and basis
         * - creates MD5 hash from basisStr + latticeStr + salt
         * @param salt Salt for hashing, empty string by default.
         * @param isScaled Whether to scale the lattice parameter 'a' to 1.
         */
        calculateHash(salt = "", isScaled = false, bypassNonPeriodicCheck = false) {
            let message;
            if (!this.isNonPeriodic || bypassNonPeriodicCheck) {
                message =
                    this.Basis.hashString + "#" + this.Lattice.getHashString(isScaled) + "#" + salt;
            }
            else {
                message = this.getInchiStringForHash();
            }
            return crypto_js_1.default.MD5(message).toString();
        }
        set hash(hash) {
            this.setProp("hash", hash);
        }
        get hash() {
            return this.prop("hash");
        }
        /**
         * Calculates hash from basis and lattice as above + scales lattice properties to make lattice.a = 1
         */
        get scaledHash() {
            return this.calculateHash("", true);
        }
        /**
         * Converts basis to crystal/fractional coordinates.
         */
        toCrystal() {
            const basis = this.Basis;
            basis.toCrystal();
            this.setProp("basis", basis.toJSON());
        }
        /**
         * Converts current material's basis coordinates to cartesian.
         * No changes if coordinates already cartesian.
         */
        toCartesian() {
            const basis = this.Basis;
            basis.toCartesian();
            this.setProp("basis", basis.toJSON());
        }
        /**
         * Returns material's basis in XYZ format.
         */
        getBasisAsXyz(fractional = false) {
            return parsers_1.default.xyz.fromMaterial(this.toJSON(), fractional);
        }
        /**
         * Returns material in Quantum Espresso output format:
         * ```
         *    CELL_PARAMETERS (angstroms)
         *    -0.543131284  -0.000000000   0.543131284
         *    -0.000000000   0.543131284   0.543131284
         *    -0.543131284   0.543131284   0.000000000
         *
         *    ATOMIC_POSITIONS (crystal)
         *    Si       0.000000000   0.000000000  -0.000000000
         *    Si       0.250000000   0.250000000   0.250000000
         * ```
         */
        getAsQEFormat() {
            return parsers_1.default.espresso.toEspressoFormat(this.toJSON());
        }
        /**
         * Returns material in POSCAR format. Pass `true` to ignore original poscar source and re-serialize.
         */
        getAsPOSCAR(ignoreOriginal = false, omitConstraints = false) {
            var _a;
            // By default return original source if exists
            if (((_a = this.src) === null || _a === void 0 ? void 0 : _a.extension) === "poscar" && !ignoreOriginal) {
                return this.src.text;
            }
            return parsers_1.default.poscar.toPoscar(this.toJSON(), omitConstraints);
        }
        /**
         * Returns a copy of the material with conventional cell constructed instead of primitive.
         */
        getACopyWithConventionalCell() {
            const material = this.clone();
            // if conventional and primitive cells are the same => return a copy.
            if ((0, conventional_cell_1.isConventionalCellSameAsPrimitiveForLatticeType)(this.Lattice.type))
                return material;
            const conventionalSupercellMatrix = conventional_cell_1.PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS[this.Lattice.type];
            const conventionalLatticeType = conventional_cell_1.PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES[this.Lattice.type];
            const config = supercell_1.default.generateConfig(material, conventionalSupercellMatrix);
            config.lattice.type = conventionalLatticeType;
            config.name = `${material.name} - conventional cell`;
            // @ts-ignore
            return new this.constructor(config);
        }
        /**
         * @summary a series of checks for the material and returns an array of results in ConsistencyChecks format.
         * @returns Array of checks results
         */
        getConsistencyChecks() {
            const basisChecks = this.getBasisConsistencyChecks();
            // any other Material checks can be added here
            return basisChecks;
        }
        /**
         * @summary a series of checks for the material's basis and returns an array of results in ConsistencyChecks format.
         * @returns Array of checks results
         */
        getBasisConsistencyChecks() {
            const checks = [];
            const limit = 1000;
            const basis = this.Basis;
            if (this.Basis.elements.length < limit) {
                const overlappingAtomsGroups = basis.getOverlappingAtoms();
                overlappingAtomsGroups.forEach(({ id1, id2, element1, element2 }) => {
                    checks.push({
                        key: `basis.coordinates.${id1}`,
                        name: "atomsOverlap",
                        severity: "warning",
                        message: `Atom ${element1} is too close to ${element2} at position ${id2 + 1}`,
                    }, {
                        key: `basis.coordinates.${id2}`,
                        name: "atomsOverlap",
                        severity: "warning",
                        message: `Atom ${element2} is too close to ${element1} at position ${id1 + 1}`,
                    });
                });
            }
            return checks;
        }
        /**
         * @summary Returns a config to create a material from a CIF or POSCAR file.
         */
        static constructMaterialFileSource(fileName, fileContent, fileExtension) {
            return {
                extension: fileExtension,
                filename: fileName,
                text: fileContent,
                hash: crypto_js_1.default.MD5(fileContent).toString(),
            };
        }
    }
    return MadeMaterial;
}
exports.MaterialMixin = MaterialMixin;
exports.Material = MaterialMixin(entity_1.HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity);
