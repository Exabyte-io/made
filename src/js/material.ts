import { type DefaultableInMemoryEntity, InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import { defaultableEntityMixin } from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import {
    type HasConsistencyChecks,
    hasConsistencyChecksMixin,
} from "@mat3ra/code/dist/js/entity/mixins/HasConsistencyChecksMixin";
import {
    type HasMetadata,
    hasMetadataMixin,
} from "@mat3ra/code/dist/js/entity/mixins/HasMetadataMixin";
import type { AnyObject } from "@mat3ra/esse/dist/js/esse/types";
import type {
    AtomicConstraintsSchema,
    ConsistencyCheck,
    DerivedPropertiesSchema,
    FileSourceSchema,
    InChIRepresentationSchema,
    LatticeSchema,
    MaterialSchema,
} from "@mat3ra/esse/dist/js/types";
import CryptoJS from "crypto-js";

import type { BasisConfig } from "./basis/basis";
import { type ConstrainedBasisConfig, ConstrainedBasis } from "./basis/constrained_basis";
import {
    isConventionalCellSameAsPrimitiveForLatticeType,
    PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES,
    PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS,
} from "./cell/conventional_cell";
import { Constraint } from "./constraints/constraints";
import { Lattice } from "./lattice/lattice";
import parsers from "./parsers/parsers";
import supercellTools from "./tools/supercell";

interface Material
    extends HasConsistencyChecks,
        DefaultableInMemoryEntity,
        HasMetadata<MaterialSchema["metadata"]> {}

// TODO: remove in-line type creation
type MaterialSchemaWithConsistencyChecksAsString = Omit<MaterialSchema, "consistencyChecks"> & {
    consistencyChecks?: ConsistencyCheck[];
};

type OptionallyConstrainedBasisConfig = BasisConfig &
    Partial<Pick<ConstrainedBasisConfig, "constraints">>;

type Schema = MaterialSchemaWithConsistencyChecksAsString;

export const defaultMaterialConfig: MaterialSchema = {
    name: "Silicon FCC",
    basis: {
        elements: [
            {
                id: 0,
                value: "Si",
            },
            {
                id: 1,
                value: "Si",
            },
        ],
        coordinates: [
            {
                id: 0,
                value: [0.0, 0.0, 0.0],
            },
            {
                id: 1,
                value: [0.25, 0.25, 0.25],
            },
        ],
        units: "crystal",
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
            length: "angstrom",
            angle: "degree",
        },
    },
};

class Material extends InMemoryEntity implements Schema {
    declare static createDefault: () => Material;

    static get defaultConfig() {
        return defaultMaterialConfig;
    }

    static constructMaterialFileSource(
        fileName: string,
        fileContent: string,
        fileExtension: string,
    ): FileSourceSchema {
        return {
            extension: fileExtension,
            filename: fileName,
            text: fileContent,
            hash: CryptoJS.MD5(fileContent).toString(),
        };
    }

    get name() {
        return this.prop("name") ?? this.formula; // if name is not set, use formula, but keep empty string
    }

    set name(name: string) {
        this.setProp("name", name);
    }

    get src() {
        return this.prop("src");
    }

    set src(src: FileSourceSchema | undefined) {
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
    getDerivedPropertyByName(name: string) {
        return this.getDerivedProperties().find((x) => x.name === name);
    }

    /**
     * @summary Returns the derived properties array for a material.
     */
    getDerivedProperties(): DerivedPropertiesSchema {
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
    setBasis(textOrObject: string | BasisConfig, format?: string, unitz?: string) {
        let basis: BasisConfig | undefined;
        if (typeof textOrObject === "string" && format === "xyz") {
            basis = parsers.xyz.toBasisConfig(textOrObject, unitz);
        } else {
            basis = textOrObject as BasisConfig;
        }
        this.setProp("basis", basis);
        this.unsetFileProps();
        this.updateFormula();
    }

    setBasisConstraints(constraints: Constraint[]) {
        const basisWithConstraints = {
            ...this.basis,
            constraints: constraints.map((c) => c.toJSON()),
        };
        this.setBasis(basisWithConstraints);
    }

    setBasisConstraintsFromArrayOfObjects(constraints: AtomicConstraintsSchema) {
        const constraintsInstances = constraints.map((c) => {
            return Constraint.fromValueAndId(c.value, c.id);
        });
        this.setBasisConstraints(constraintsInstances);
    }

    get basis() {
        return this.requiredProp<OptionallyConstrainedBasisConfig>("basis");
    }

    // returns the instance of {ConstrainedBasis} class
    get Basis() {
        const basisData = this.basis;

        return new ConstrainedBasis({
            ...basisData,
            cell: this.Lattice.vectors,
            constraints: basisData.constraints || [],
        });
    }

    /**
     * High-level access to unique elements from material instead of basis.
     */
    get uniqueElements() {
        return this.Basis.uniqueElements;
    }

    get lattice(): LatticeSchema {
        return this.prop("lattice") as LatticeSchema;
    }

    set lattice(config: LatticeSchema) {
        const originalIsInCrystalUnits = this.Basis.isInCrystalUnits;
        const basis = this.Basis;
        basis.toCartesian();

        const newLattice = new Lattice(config);
        basis.cell = newLattice.vectors;
        if (originalIsInCrystalUnits) basis.toCrystal();

        // Preserve all properties from the original basis to ensure constraints are included
        const newBasisConfig = basis.toJSON();
        this.setProp("basis", newBasisConfig);
        this.setProp("lattice", config);

        this.unsetFileProps();
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
            return (inchi as InChIRepresentationSchema).value;
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
        return this.prop("hash") as string;
    }

    /**
     * Calculates hash from basis and lattice as above + scales lattice properties to make lattice.a = 1
     */
    get scaledHash(): string {
        return this.calculateHash("", true);
    }

    get external() {
        return this.prop<MaterialSchema["external"]>("external");
    }

    set external(external: MaterialSchema["external"]) {
        this.setProp("external", external);
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
        // By default return original source if exists
        if (this.src?.extension === "poscar" && !ignoreOriginal) {
            return this.src.text;
        }
        return parsers.poscar.toPoscar(this.toJSON(), omitConstraints);
    }

    /**
     * Returns a copy of the material with conventional cell constructed instead of primitive.
     */
    getACopyWithConventionalCell(): this {
        const material = this.clone();

        // if conventional and primitive cells are the same => return a copy.
        if (isConventionalCellSameAsPrimitiveForLatticeType(this.Lattice.type)) {
            return material;
        }

        const conventionalSupercellMatrix =
            PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS[this.Lattice.type];
        const conventionalLatticeType =
            PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES[this.Lattice.type];
        const config = supercellTools.generateConfig(this, conventionalSupercellMatrix);

        config.lattice.type = conventionalLatticeType;
        config.name = `${this.name} - conventional cell`;

        // @ts-expect-error
        return new this.constructor(config) as this;
    }

    /**
     * @summary a series of checks for the material and returns an array of results in ConsistencyChecks format.
     * @returns Array of checks results
     */
    getConsistencyChecks(): ConsistencyCheck[] {
        const basisChecks = this.getBasisConsistencyChecks();

        // any other Material checks can be added here

        return basisChecks;
    }

    /**
     * @summary a series of checks for the material's basis and returns an array of results in ConsistencyChecks format.
     * @returns Array of checks results
     */
    getBasisConsistencyChecks(): ConsistencyCheck[] {
        const checks: ConsistencyCheck[] = [];
        const limit = 1000;
        const basis = this.Basis;

        if (this.Basis.elements.length < limit) {
            const overlappingAtomsGroups = basis.getOverlappingAtoms();
            overlappingAtomsGroups.forEach(({ id1, id2, element1, element2 }) => {
                checks.push(
                    {
                        key: `basis.coordinates.${id1}`,
                        name: "atomsOverlap",
                        severity: "warning",
                        message: `Atom ${element1} is too close to ${element2} at position ${
                            id2 + 1
                        }`,
                    },
                    {
                        key: `basis.coordinates.${id2}`,
                        name: "atomsOverlap",
                        severity: "warning",
                        message: `Atom ${element2} is too close to ${element1} at position ${
                            id1 + 1
                        }`,
                    },
                );
            });
        }

        return checks;
    }

    toJSON(): MaterialSchema & AnyObject {
        // validation intoJSON() method will fail if name is not set in _json
        this.name = this.name || this.formula;

        return {
            ...(super.toJSON() as MaterialSchema & AnyObject),
            lattice: this.Lattice.toJSON(),
            basis: this.Basis.toJSON(),
            isNonPeriodic: this.isNonPeriodic,
        };
    }
}

defaultableEntityMixin(Material);
hasConsistencyChecksMixin(Material.prototype);
hasMetadataMixin(Material.prototype);

export { Material };
