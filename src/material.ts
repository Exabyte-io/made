import { HasMetadataNamedDefaultableInMemoryEntity } from "@exabyte-io/code.js/dist/entity";
import CryptoJS from "crypto-js";

import {
    MaterialSchema,
} from "@exabyte-io/code.js/src/types"


import { ConstrainedBasis, ConstrainedBasisJSON } from "./basis/constrained_basis";
import {
    isConventionalCellSameAsPrimitiveForLatticeType,
    PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES,
    PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS,
} from "./cell/conventional_cell";
import { ATOMIC_COORD_UNITS, units } from "./constants";
import { Lattice, LatticeJSON } from "./lattice/lattice";
import { BravaisConfigProps } from "./lattice/lattice_vectors";
import { LATTICE_TYPE } from "./lattice/types";
import parsers from "./parsers/parsers";
import supercellTools from "./tools/supercell";
import { BasisConfig } from "./parsers/xyz";
import { Constraint } from "./constraints/constraints";
import { AnyObject } from "@exabyte-io/code.js/dist/entity/in_memory";

export const defaultMaterialConfig = {
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
        units: ATOMIC_COORD_UNITS.crystal,
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
            angle: units.degree,
        },
    },
};

interface Property {
    name: string;
    value: string;
}

export interface MaterialJSON extends AnyObject {
    lattice:LatticeJSON;
    basis: ConstrainedBasisJSON;
    name: string;
    isNonPeriodic: boolean;
}

interface MaterialSchemaJSON extends MaterialSchema, AnyObject {

}


export abstract class Material extends HasMetadataNamedDefaultableInMemoryEntity {
    abstract src: {
        extension: string;
        text: string;
    };

    _json: MaterialSchemaJSON;

    constructor(config: MaterialSchemaJSON) {
        super(config);
        this._json = {...config};
        this.name = super.name || this.formula
    }

    toJSON(): MaterialJSON {
        return {
            lattice: this.Lattice.toJSON(),
            basis: this.Basis.toJSON(),
            name: this.name,
            isNonPeriodic: this.isNonPeriodic,
        };
    }

    static get defaultConfig() {
        return defaultMaterialConfig;
    }

    updateFormula() {
        this.setProp("formula", this.Basis.formula);
        this.setProp("unitCellFormula", this.Basis.unitCellFormula);
    }

    /**
     * Gets Bolean value for whether or not a material is non-periodic vs periodic.
     * False = periodic, True = non-periodic
     */
    get isNonPeriodic(): boolean {
        return this.prop("isNonPeriodic", false);
    }

    /**
     * @summary Sets the value of isNonPeriodic based on Boolean value passed as an argument.
     */
    set isNonPeriodic(bool: boolean) {
        this.setProp("isNonPeriodic", bool);
    }

    /**
     * @summary Returns the specific derived property (as specified by name) for a material.
     */
    getDerivedPropertyByName(name: string): Property| undefined {
        return this.getDerivedProperties().find((x) => x.name === name);
    }

    /**
     * @summary Returns the derived properties array for a material.
     */
    getDerivedProperties(): Property[] {
        return this.prop("derivedProperties", []);
    }

    /**
     * Gets material's formula
     */
    get formula(): string {
        return this.prop("formula") || this.Basis.formula;
    }

    get unitCellFormula(): string {
        return this.prop("unitCellFormula") || this.Basis.unitCellFormula;
    }

    /**
     * @param textOrObject Basis text or JSON object.
     * @param format Format (xyz, etc.)
     * @param unitz crystal/cartesian
     */
    setBasis(textOrObject: string | BasisConfig, format?: string, unitz?: string) {
        let basis: BasisConfig | undefined;
        switch (format) {
            case "xyz":
                basis = parsers.xyz.toBasisConfig(textOrObject as string, unitz);
                break;
            default:
                basis = textOrObject as BasisConfig;
        }
        this.setProp("basis", basis);
        this.updateFormula();
    }

    setBasisConstraints(constraints: Constraint[]) {
        this.setBasis({ ...this.basis, constraints });
    }

    get basis(): BasisConfig {
        return this.prop<BasisConfig>("basis", undefined);
    }

    // returns the instance of {ConstrainedBasis} class
    get Basis() {
        return new ConstrainedBasis({
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

    get lattice(): BravaisConfigProps | undefined {
        return this.prop("lattice", undefined);
    }

    set lattice(config: BravaisConfigProps | undefined) {
        this.setProp("lattice", config);
    }

    get Lattice(): Lattice {
        return new Lattice(this.lattice);
    }

    /**
     * Returns the inchi string from the derivedProperties for a non-periodic material, or throws an error if the
     *  inchi cannot be found.
     *  @returns {String}
     */
    getInchiStringForHash(): string {
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
    calculateHash(salt = "", isScaled = false, bypassNonPeriodicCheck = false): string {
        let message;
        if (!this.isNonPeriodic || bypassNonPeriodicCheck) {
            message =
                this.Basis.hashString + "#" + this.Lattice.getHashString(isScaled) + "#" + salt;
        } else {
            message = this.getInchiStringForHash();
        }
        return CryptoJS.MD5(message).toString();
    }

    set hash(hash: string) {
        this.setProp("hash", hash);
    }

    get hash(): string {
        return this.prop("hash");
    }

    /**
     * Calculates hash from basis and lattice as above + scales lattice properties to make lattice.a = 1
     */
    get scaledHash(): string {
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
    getBasisAsXyz(fractional = false): string {
        return parsers.xyz.fromMaterial(this.toJSON(), fractional);
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
    getAsQEFormat(): string {
        return parsers.espresso.toEspressoFormat(this.toJSON());
    }

    /**
     * Returns material in POSCAR format. Pass `true` to ignore original poscar source and re-serialize.
     */
    getAsPOSCAR(ignoreOriginal = false, omitConstraints = false): string {
        const { src } = this;
        // By default return original source if exists
        if (src && src.extension === "poscar" && !ignoreOriginal) {
            return this.src.text;
        }
        return parsers.poscar.toPoscar(this.toJSON(), omitConstraints);
    }

    /**
     * Returns a copy of the material with conventional cell constructed instead of primitive.
     */
    getACopyWithConventionalCell() {
        const material = this.clone<Material>();

        // if conventional and primitive cells are the same => return a copy.
        if (isConventionalCellSameAsPrimitiveForLatticeType(this.Lattice.type)) return material;

        const conventionalSupercellMatrix =
            PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS[this.Lattice.type];
        const conventionalLatticeType =
            PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES[this.Lattice.type];
        const config = supercellTools.generateConfig(material, conventionalSupercellMatrix);

        config.lattice.type = conventionalLatticeType;
        config.name = `${material.name} - conventional cell`;

        // @ts-ignore
        return new this.constructor(config);
    }
}
