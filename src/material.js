import lodash from "lodash";
import {mix} from "mixwith";
import CryptoJS from "crypto-js";

import parsers from "./parsers/parsers";
import {Lattice} from "./lattice/lattice";
import {LATTICE_TYPE} from "./lattice/types";
import {ATOMIC_COORD_UNITS, units} from "./constants";
import {ConstrainedBasis} from "./basis/constrained_basis";

export const defaultMaterialConfig = {
    name: 'Silicon FCC',
    basis: {
        elements: [
            {
                id: 1,
                value: 'Si'
            },
            {
                id: 2,
                value: 'Si'
            }
        ],
        coordinates: [
            {
                id: 1,
                value: [
                    0.00,
                    0.00,
                    0.00
                ]
            },
            {
                id: 2,
                value: [
                    0.25,
                    0.25,
                    0.25
                ]
            }
        ],
        units: ATOMIC_COORD_UNITS.crystal
    },
    lattice: {
        // Primitive cell for Diamond FCC Silicon at ambient conditions
        type: LATTICE_TYPE.FCC,
        a: 3.867,
        b: 3.867,
        c: 3.867,
        alpha: 60,
        beta: 60,
        gamma: 60,
        units: {
            length: units.angstrom,
            angle: units.degree
        }
    }
};

export class Material {
    constructor(config) {
        this._json = lodash.cloneDeep(config || {});
    }

    prop(name, defaultValue) {
        // `lodash.get` gets `null` when the value is `null`, but we still want a default value in this case, hence `||`
        return lodash.get(this._json, name, defaultValue) || defaultValue;
    }

    setProp(name, value) {
        this._json[name] = value;
    }

    static get defaultConfig() {
        return defaultMaterialConfig;
    }

    _deepClone(object) {
        return lodash.cloneDeep(object);
    }

    updateFormula() {
        this.setProp('formula', this.Basis.formula);
        this.setProp('unitCellFormula', this.Basis.unitCellFormula);
    }

    /**
     * @summary Gets material's formula
     */
    get formula() {
        return this.prop('formula') || this.Basis.formula;
    }

    get unitCellFormula() {
        return this.prop('unitCellFormula') || this.Basis.unitCellFormula;
    }

    get name() {return this.prop('name') || this.formula}

    set name(name) {this.setProp('name', name)}

    /**
     * @param textOrObject {String} Basis text or JSON object.
     * @param format {String} Format (xyz, etc.)
     * @param units {String} crystal/cartesian
     */
    setBasis(textOrObject, format, units) {
        let basis;
        switch (format) {
            case 'xyz':
                basis = parsers.xyz.toBasisConfig(textOrObject, units);
                break;
            default:
                basis = textOrObject;

        }
        this.setProp('basis', basis);
        this.updateFormula();
    }

    setBasisConstraints(constraints) {
        this.setBasis({
            ...this.basis,
            constraints,
        })
    }

    clone(extraContext) {
        return new this.constructor(Object.assign({}, this.toJSON(), extraContext));
    }

    get basis() {
        return this.prop('basis', undefined, true);
    }

    // returns the instance of {ConstrainedBasis} class
    get Basis() {
        return new ConstrainedBasis({
            ...this.basis,
            cell: this.Lattice.vectorArrays
        });
    }

    get lattice() {
        return this.prop('lattice', undefined, true);
    }

    set lattice(config) {
        return this.setProp('lattice', config);
    }

    // returns the instance of {Lattice} class
    get Lattice() {
        return new Lattice(this.lattice);
    }

    /**
     * @summary Calculates hash from basis and lattice. Algorithm expects the following:
     * - asserts lattice units to be angstrom
     * - asserts basis units to be crystal
     * - asserts basis coordinates and lattice measurements are rounded to hash precision
     * - forms strings for lattice and basis
     * - creates MD5 hash from basisStr + latticeStr + salt
     * @param salt {String} Salt for hashing, empty string by default.
     * @param isScaled {Boolean} Whether to scale the lattice parameter 'a' to 1.
     */
    calculateHash(salt = '', isScaled = false) {
        const message = this.Basis.hashString + "#" + this.Lattice.getHashString(isScaled) + "#" + salt;
        return CryptoJS.MD5(message).toString();
    }

    get hash() {
        return this.calculateHash();
    }

    /**
     * @summary Calculates hash from basis and lattice as above + scales lattice properties to make lattice.a = 1
     */
    get scaledHash() {
        return this.calculateHash('', true);
    }

    /**
     * @summary Converts basis to crystal/fractional coordinates.
     */
    toCrystal() {
        const basis = this.Basis;
        basis.toCrystal();
        this.setProp('basis', basis.toJSON())
    }

    /**
     * @summary Converts current material's basis coordinates to cartesian.
     * No changes if coordinates already cartesian.
     */
    toCartesian() {
        const basis = this.Basis;
        basis.toCartesian();
        this.setProp('basis', basis.toJSON())
    }

    /**
     * @summary Returns material's basis in XYZ format.
     * @return {String}
     */
    getBasisAsXyz(fractional = false) {
        return parsers.xyz.fromMaterial(this.toJSON(), fractional);
    }

    /**
     * @summary Returns material in Quantum Espresso output format:
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
     * @return {String}
     */
    getAsQEFormat() {
        return parsers.espresso.toEspressoFormat(this.toJSON());
    }

    /**
     * @summary Returns material in POSCAR format. Pass `true` to ignore original poscar source and re-serialize.
     */
    getAsPOSCAR(ignoreOriginal = false, omitConstraints = false) {
        const src = this.src;
        // By default return original source if exists
        if ((src && src.extension === 'poscar') && !ignoreOriginal) {
            return this.src.text;
        }
        return parsers.poscar.toPoscar(this.toJSON(), omitConstraints);
    }

    toJSON() {
        return {
            lattice: this.Lattice.toJSON(),
            basis: this.Lattice.toJSON(),
            name: this.name || this.formula
        };
    }

}
