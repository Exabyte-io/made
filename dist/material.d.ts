import { HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity } from "@exabyte-io/code.js/dist/entity";
import { AnyObject } from "@exabyte-io/code.js/dist/entity/in_memory";
import { ConsistencyCheck, DerivedPropertiesSchema, FileSourceSchema, MaterialSchema } from "@mat3ra/esse/lib/js/types";
import { ConstrainedBasis } from "./basis/constrained_basis";
import { Constraint } from "./constraints/constraints";
import { Lattice } from "./lattice/lattice";
import { BravaisConfigProps } from "./lattice/lattice_vectors";
import { BasisConfig } from "./parsers/xyz";
import { MaterialJSON } from "./types";
export declare const defaultMaterialConfig: {
    name: string;
    basis: {
        elements: {
            id: number;
            value: string;
        }[];
        coordinates: {
            id: number;
            value: number[];
        }[];
        units: string;
    };
    lattice: {
        type: string;
        a: number;
        b: number;
        c: number;
        alpha: number;
        beta: number;
        gamma: number;
        units: {
            length: string;
            angle: string;
        };
    };
};
export interface MaterialSchemaJSON extends MaterialSchema, AnyObject {
}
type MaterialBaseEntity = InstanceType<typeof HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity>;
export type MaterialBaseEntityConstructor<T extends MaterialBaseEntity = MaterialBaseEntity> = new (...args: any[]) => T;
export declare function MaterialMixin<T extends MaterialBaseEntityConstructor = MaterialBaseEntityConstructor>(superclass: T): {
    new (...config: any[]): {
        _json: MaterialSchemaJSON;
        toJSON(): MaterialJSON;
        src: FileSourceSchema;
        updateFormula(): void;
        /**
         * Gets Bolean value for whether or not a material is non-periodic vs periodic.
         * False = periodic, True = non-periodic
         */
        isNonPeriodic: boolean;
        /**
         * @summary Returns the specific derived property (as specified by name) for a material.
         */
        getDerivedPropertyByName(name: string): {
            name?: "volume" | undefined;
            units?: "angstrom^3" | undefined;
            value: number;
        } | {
            name?: "density" | undefined;
            units?: "g/cm^3" | undefined;
            value: number;
        } | {
            pointGroupSymbol?: string | undefined;
            spaceGroupSymbol?: string | undefined;
            tolerance?: {
                units?: "angstrom" | undefined;
                value: number;
            } | undefined;
            name?: "symmetry" | undefined;
        } | {
            name?: "elemental_ratio" | undefined;
            value: number;
            element?: string | undefined;
        } | {
            name?: "p-norm" | undefined;
            degree?: number | undefined;
            value: number;
        } | {
            name?: "inchi" | undefined;
            value: string;
        } | {
            name?: "inchi_key" | undefined;
            value: string;
        } | undefined;
        /**
         * @summary Returns the derived properties array for a material.
         */
        getDerivedProperties(): DerivedPropertiesSchema;
        /**
         * Gets material's formula
         */
        readonly formula: string;
        readonly unitCellFormula: string;
        /**
         * @param textOrObject Basis text or JSON object.
         * @param format Format (xyz, etc.)
         * @param unitz crystal/cartesian
         */
        setBasis(textOrObject: string | BasisConfig, format?: string, unitz?: string): void;
        setBasisConstraints(constraints: Constraint[]): void;
        readonly basis: BasisConfig;
        readonly Basis: ConstrainedBasis;
        /**
         * High-level access to unique elements from material instead of basis.
         */
        readonly uniqueElements: string[];
        lattice: BravaisConfigProps | undefined;
        readonly Lattice: Lattice;
        /**
         * Returns the inchi string from the derivedProperties for a non-periodic material, or throws an error if the
         *  inchi cannot be found.
         *  @returns {String}
         */
        getInchiStringForHash(): string;
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
        calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
        hash: string;
        /**
         * Calculates hash from basis and lattice as above + scales lattice properties to make lattice.a = 1
         */
        readonly scaledHash: string;
        /**
         * Converts basis to crystal/fractional coordinates.
         */
        toCrystal(): void;
        /**
         * Converts current material's basis coordinates to cartesian.
         * No changes if coordinates already cartesian.
         */
        toCartesian(): void;
        /**
         * Returns material's basis in XYZ format.
         */
        getBasisAsXyz(fractional?: boolean): string;
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
        getAsQEFormat(): string;
        /**
         * Returns material in POSCAR format. Pass `true` to ignore original poscar source and re-serialize.
         */
        getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
        /**
         * Returns a copy of the material with conventional cell constructed instead of primitive.
         */
        getACopyWithConventionalCell(): any;
        /**
         * @summary a series of checks for the material and returns an array of results in ConsistencyChecks format.
         * @returns Array of checks results
         */
        getConsistencyChecks(): ConsistencyCheck[];
        /**
         * @summary a series of checks for the material's basis and returns an array of results in ConsistencyChecks format.
         * @returns Array of checks results
         */
        getBasisConsistencyChecks(): ConsistencyCheck[];
        consistencyChecks: object[];
        addConsistencyChecks(array: object[]): void;
        prop: {
            <T_1 = undefined>(name: string, defaultValue: T_1): T_1;
            <T_1 = undefined>(name: string): T_1 | undefined;
        } & {
            <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
            <T_1_1 = undefined>(name: string): T_1_1 | undefined;
        } & {
            <T_3 = undefined>(name: string, defaultValue: T_3): T_3;
            <T_1_2 = undefined>(name: string): T_1_2 | undefined;
        } & {
            <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
            <T_1_3 = undefined>(name: string): T_1_3 | undefined;
        } & {
            <T_5 = undefined>(name: string, defaultValue: T_5): T_5;
            <T_6 = undefined>(name: string): T_6 | undefined;
        };
        setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
        unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
        toJSONSafe: ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject);
        toJSONQuick: ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject);
        clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
        validate: (() => void) & (() => void) & (() => void) & (() => void) & (() => void);
        clean: ((config: AnyObject) => AnyObject) & ((config: AnyObject) => AnyObject) & ((config: AnyObject) => AnyObject) & ((config: AnyObject) => AnyObject) & ((config: AnyObject) => AnyObject);
        isValid: (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean);
        id: string;
        readonly cls: string;
        getClsName: (() => string) & (() => string) & (() => string) & (() => string) & (() => string);
        readonly slug: string;
        readonly isSystemEntity: boolean;
        getAsEntityReference: ((byIdOnly?: boolean | undefined) => import("@mat3ra/esse/lib/js/types").EntityReferenceSchema) & ((byIdOnly?: boolean | undefined) => import("@mat3ra/esse/lib/js/types").EntityReferenceSchema) & ((byIdOnly?: boolean | undefined) => import("@mat3ra/esse/lib/js/types").EntityReferenceSchema) & ((byIdOnly?: boolean | undefined) => import("@mat3ra/esse/lib/js/types").EntityReferenceSchema) & ((byIdOnly?: boolean | undefined) => import("@mat3ra/esse/lib/js/types").EntityReferenceSchema);
        getEntityByName: ((entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string) => import("@exabyte-io/code.js/dist/entity").InMemoryEntity) & ((entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string) => import("@exabyte-io/code.js/dist/entity").InMemoryEntity) & ((entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string) => import("@exabyte-io/code.js/dist/entity").InMemoryEntity) & ((entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string) => import("@exabyte-io/code.js/dist/entity").InMemoryEntity) & ((entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string) => import("@exabyte-io/code.js/dist/entity").InMemoryEntity);
        metadata: object;
        updateMetadata(object: object): void;
        name: string;
        setName(name: string): void;
        readonly isDefault: boolean;
    };
    readonly jsonSchema: {
        $id: string;
        $schema: string;
        title: string;
        type: string;
        required: string[];
        properties: {
            formula: {
                description: string;
                type: string;
            };
            unitCellFormula: {
                description: string;
                type: string;
            };
            basis: {
                $schema: string;
                title: string;
                type: string;
                required: string[];
                properties: {
                    elements: {
                        type: string;
                        items: {
                            $schema: string;
                            title: string;
                            description: string;
                            type: string;
                            required: string[];
                            properties: {
                                id: {
                                    type: string;
                                };
                                value: {
                                    type: string;
                                };
                                occurrence: {
                                    description: string;
                                    type: string;
                                };
                                oxidationState: {
                                    type: string;
                                };
                            };
                        };
                    };
                    labels: {
                        description: string;
                        type: string;
                        items: {
                            properties: {
                                id: {
                                    type: string;
                                };
                                value: {
                                    type: string;
                                };
                            };
                        };
                    };
                    coordinates: {
                        type: string;
                        items: {
                            $schema: string;
                            title: string;
                            description: string;
                            type: string;
                            properties: {
                                id: {
                                    type: string;
                                };
                                value: {
                                    $schema: string;
                                    title: string;
                                    oneOf: {
                                        $schema: string;
                                        title: string;
                                        type: string;
                                        items: {
                                            type: string;
                                        };
                                        minItems: number;
                                        maxItems: number;
                                    }[];
                                };
                            };
                        };
                    };
                    name: {
                        type: string;
                    };
                    units: {
                        type: string;
                    };
                    bonds: {
                        $schema: string;
                        title: string;
                        type: string;
                        uniqueItems: boolean;
                        items: {
                            type: string;
                            properties: {
                                atomPair: {
                                    description: string;
                                    type: string;
                                    minItems: number;
                                    maxItems: number;
                                    $schema: string;
                                    title: string;
                                    items: {
                                        type: string;
                                        properties: {
                                            id: {
                                                description: string;
                                                type: string;
                                            };
                                        };
                                    };
                                };
                                bondType: {
                                    type: string;
                                    enum: string[];
                                };
                            }; /**
                             * Returns the inchi string from the derivedProperties for a non-periodic material, or throws an error if the
                             *  inchi cannot be found.
                             *  @returns {String}
                             */
                        };
                    };
                };
            };
            lattice: {
                $schema: string;
                title: string;
                type: string;
                required: string[];
                properties: {
                    name: {
                        enum: string[];
                    };
                    vectors: {
                        $schema: string;
                        title: string;
                        type: string;
                        required: string[];
                        properties: {
                            alat: {
                                description: string;
                                type: string;
                                default: number;
                            };
                            units: {
                                enum: string[];
                            };
                            a: {
                                $schema: string;
                                title: string;
                                type: string;
                                minItems: number;
                                maxItems: number;
                                items: {
                                    type: string;
                                };
                            };
                            b: {
                                $schema: string;
                                title: string;
                                type: string;
                                minItems: number;
                                maxItems: number;
                                items: {
                                    type: string;
                                };
                            };
                            c: {
                                $schema: string;
                                title: string;
                                type: string;
                                minItems: number;
                                maxItems: number;
                                items: {
                                    type: string;
                                };
                            };
                        };
                    };
                    type: {
                        $schema: string;
                        title: string;
                        type: string;
                        enum: string[];
                    };
                    units: {
                        type: string;
                        properties: {
                            length: {
                                type: string;
                                enum: string[];
                            };
                            angle: {
                                type: string;
                                enum: string[];
                            };
                        };
                    };
                    a: {
                        description: string;
                        type: string;
                    };
                    b: {
                        description: string;
                        type: string;
                    };
                    c: {
                        description: string;
                        type: string;
                    };
                    alpha: {
                        description: string;
                        type: string;
                    };
                    beta: {
                        description: string;
                        type: string;
                    };
                    gamma: {
                        description: string;
                        type: string;
                    };
                };
            };
            derivedProperties: {
                $schema: string;
                title: string;
                type: string;
                items: {
                    oneOf: ({
                        $schema: string;
                        title: string;
                        type: string;
                        required: string[];
                        properties: {
                            name: {
                                enum: string[];
                            };
                            units: {
                                enum: string[];
                            };
                            value: {
                                type: string;
                                minimum?: undefined;
                                maximum?: undefined;
                            };
                            pointGroupSymbol?: undefined;
                            spaceGroupSymbol?: undefined;
                            tolerance?: undefined;
                            element?: undefined;
                            degree?: undefined;
                        };
                        description?: undefined;
                    } | {
                        $schema: string;
                        title: string;
                        type: string;
                        properties: {
                            pointGroupSymbol: {
                                description: string;
                                type: string;
                            };
                            spaceGroupSymbol: {
                                description: string;
                                type: string;
                            };
                            tolerance: {
                                type: string;
                                description: string;
                                $schema: string;
                                title: string;
                                required: string[];
                                properties: {
                                    units: {
                                        enum: string[];
                                    };
                                    value: {
                                        type: string;
                                    };
                                };
                            };
                            name: {
                                enum: string[];
                            };
                            units?: undefined;
                            value?: undefined;
                            element?: undefined;
                            degree?: undefined;
                        };
                        required?: undefined;
                        description?: undefined;
                    } | {
                        $schema: string;
                        title: string;
                        description: string;
                        type: string;
                        required: string[];
                        properties: {
                            name: {
                                enum: string[];
                            };
                            value: {
                                type: string;
                                minimum: number;
                                maximum: number;
                            };
                            element: {
                                type: string;
                                description: string;
                            };
                            units?: undefined;
                            pointGroupSymbol?: undefined;
                            spaceGroupSymbol?: undefined;
                            tolerance?: undefined;
                            degree?: undefined;
                        };
                    } | {
                        $schema: string;
                        title: string;
                        description: string;
                        type: string;
                        required: string[];
                        properties: {
                            name: {
                                enum: string[];
                            };
                            degree: {
                                type: string;
                                description: string;
                            };
                            value: {
                                type: string;
                                minimum?: undefined;
                                maximum?: undefined;
                            };
                            units?: undefined;
                            pointGroupSymbol?: undefined;
                            spaceGroupSymbol?: undefined;
                            tolerance?: undefined;
                            element?: undefined;
                        };
                    } | {
                        $schema: string;
                        title: string;
                        type: string;
                        required: string[];
                        properties: {
                            name: {
                                enum: string[];
                            };
                            value: {
                                type: string;
                                minimum?: undefined;
                                maximum?: undefined;
                            };
                            units?: undefined;
                            pointGroupSymbol?: undefined;
                            spaceGroupSymbol?: undefined;
                            tolerance?: undefined;
                            element?: undefined;
                            degree?: undefined;
                        };
                        description?: undefined;
                    })[];
                    discriminator: {
                        propertyName: string;
                    };
                    required: string[];
                };
            };
            external: {
                $schema: string;
                title: string;
                description: string;
                type: string;
                required: string[];
                properties: {
                    id: {
                        description: string;
                        oneOf: {
                            type: string;
                        }[];
                    };
                    source: {
                        description: string;
                        type: string;
                    };
                    origin: {
                        description: string;
                        type: string;
                    };
                    data: {
                        description: string;
                        type: string;
                    };
                    doi: {
                        description: string;
                        type: string;
                    };
                    url: {
                        description: string;
                        type: string;
                    };
                };
            };
            src: {
                $schema: string;
                title: string;
                description: string;
                type: string;
                required: string[];
                properties: {
                    extension: {
                        description: string;
                        type: string;
                    };
                    filename: {
                        description: string;
                        type: string;
                    };
                    text: {
                        description: string;
                        type: string;
                    };
                    hash: {
                        description: string;
                        type: string;
                    };
                };
            };
            scaledHash: {
                description: string;
                type: string;
            };
            icsdId: {
                description: string;
                type: string;
            };
            isNonPeriodic: {
                description: string;
                type: string;
            };
            _id: {
                description: string;
                type: string;
            };
            slug: {
                description: string;
                type: string;
            };
            systemName: {
                type: string;
            };
            consistencyChecks: {
                type: string;
                items: {
                    $schema: string;
                    title: string;
                    type: string;
                    description: string;
                    required: string[];
                    properties: {
                        key: {
                            type: string;
                            description: string;
                        };
                        name: {
                            enum: string[];
                            description: string;
                        };
                        severity: {
                            enum: string[];
                            description: string;
                        };
                        message: {
                            type: string;
                            description: string;
                        };
                    };
                };
            };
            schemaVersion: {
                description: string;
                type: string;
                default: string;
            };
            name: {
                description: string;
                type: string;
            };
            isDefault: {
                description: string;
                type: string;
                default: boolean;
            };
            metadata: {
                type: string;
            };
        };
    };
    readonly defaultConfig: {
        name: string;
        basis: {
            elements: {
                id: number;
                value: string;
            }[];
            coordinates: {
                id: number;
                value: number[];
            }[];
            units: string;
        };
        lattice: {
            type: string;
            a: number;
            b: number;
            c: number;
            alpha: number;
            beta: number;
            gamma: number;
            units: {
                length: string;
                angle: string;
            };
        };
    };
} & T;
export declare const Material: {
    new (...config: any[]): {
        _json: MaterialSchemaJSON;
        toJSON(): MaterialJSON;
        src: FileSourceSchema;
        updateFormula(): void;
        /**
         * Gets Bolean value for whether or not a material is non-periodic vs periodic.
         * False = periodic, True = non-periodic
         */
        isNonPeriodic: boolean;
        /**
         * @summary Returns the specific derived property (as specified by name) for a material.
         */
        getDerivedPropertyByName(name: string): {
            name?: "volume" | undefined;
            units?: "angstrom^3" | undefined;
            value: number;
        } | {
            name?: "density" | undefined;
            units?: "g/cm^3" | undefined;
            value: number;
        } | {
            pointGroupSymbol?: string | undefined;
            spaceGroupSymbol?: string | undefined;
            tolerance?: {
                units?: "angstrom" | undefined;
                value: number;
            } | undefined;
            name?: "symmetry" | undefined;
        } | {
            name?: "elemental_ratio" | undefined;
            value: number;
            element?: string | undefined;
        } | {
            name?: "p-norm" | undefined;
            degree?: number | undefined;
            value: number;
        } | {
            name?: "inchi" | undefined;
            value: string;
        } | {
            name?: "inchi_key" | undefined;
            value: string;
        } | undefined;
        /**
         * @summary Returns the derived properties array for a material.
         */
        getDerivedProperties(): DerivedPropertiesSchema;
        /**
         * Gets material's formula
         */
        readonly formula: string;
        readonly unitCellFormula: string;
        /**
         * @param textOrObject Basis text or JSON object.
         * @param format Format (xyz, etc.)
         * @param unitz crystal/cartesian
         */
        setBasis(textOrObject: string | BasisConfig, format?: string, unitz?: string): void;
        setBasisConstraints(constraints: Constraint[]): void;
        readonly basis: BasisConfig;
        readonly Basis: ConstrainedBasis;
        /**
         * High-level access to unique elements from material instead of basis.
         */
        readonly uniqueElements: string[];
        lattice: BravaisConfigProps | undefined;
        readonly Lattice: Lattice;
        /**
         * Returns the inchi string from the derivedProperties for a non-periodic material, or throws an error if the
         *  inchi cannot be found.
         *  @returns {String}
         */
        getInchiStringForHash(): string;
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
        calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
        hash: string;
        /**
         * Calculates hash from basis and lattice as above + scales lattice properties to make lattice.a = 1
         */
        readonly scaledHash: string;
        /**
         * Converts basis to crystal/fractional coordinates.
         */
        toCrystal(): void;
        /**
         * Converts current material's basis coordinates to cartesian.
         * No changes if coordinates already cartesian.
         */
        toCartesian(): void;
        /**
         * Returns material's basis in XYZ format.
         */
        getBasisAsXyz(fractional?: boolean): string;
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
        getAsQEFormat(): string;
        /**
         * Returns material in POSCAR format. Pass `true` to ignore original poscar source and re-serialize.
         */
        getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
        /**
         * Returns a copy of the material with conventional cell constructed instead of primitive.
         */
        getACopyWithConventionalCell(): any;
        /**
         * @summary a series of checks for the material and returns an array of results in ConsistencyChecks format.
         * @returns Array of checks results
         */
        getConsistencyChecks(): ConsistencyCheck[];
        /**
         * @summary a series of checks for the material's basis and returns an array of results in ConsistencyChecks format.
         * @returns Array of checks results
         */
        getBasisConsistencyChecks(): ConsistencyCheck[];
        consistencyChecks: object[];
        addConsistencyChecks(array: object[]): void;
        prop: {
            <T = undefined>(name: string, defaultValue: T): T;
            <T_1 = undefined>(name: string): T_1 | undefined;
        } & {
            <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
            <T_1_1 = undefined>(name: string): T_1_1 | undefined;
        } & {
            <T_3 = undefined>(name: string, defaultValue: T_3): T_3;
            <T_1_2 = undefined>(name: string): T_1_2 | undefined;
        } & {
            <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
            <T_1_3 = undefined>(name: string): T_1_3 | undefined;
        } & {
            <T_5 = undefined>(name: string, defaultValue: T_5): T_5;
            <T_6 = undefined>(name: string): T_6 | undefined;
        };
        setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
        unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
        toJSONSafe: ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject);
        toJSONQuick: ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject) & ((exclude?: string[] | undefined) => AnyObject);
        clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
        validate: (() => void) & (() => void) & (() => void) & (() => void) & (() => void);
        clean: ((config: AnyObject) => AnyObject) & ((config: AnyObject) => AnyObject) & ((config: AnyObject) => AnyObject) & ((config: AnyObject) => AnyObject) & ((config: AnyObject) => AnyObject);
        isValid: (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean);
        id: string;
        readonly cls: string;
        getClsName: (() => string) & (() => string) & (() => string) & (() => string) & (() => string);
        readonly slug: string;
        readonly isSystemEntity: boolean;
        getAsEntityReference: ((byIdOnly?: boolean | undefined) => import("@mat3ra/esse/lib/js/types").EntityReferenceSchema) & ((byIdOnly?: boolean | undefined) => import("@mat3ra/esse/lib/js/types").EntityReferenceSchema) & ((byIdOnly?: boolean | undefined) => import("@mat3ra/esse/lib/js/types").EntityReferenceSchema) & ((byIdOnly?: boolean | undefined) => import("@mat3ra/esse/lib/js/types").EntityReferenceSchema) & ((byIdOnly?: boolean | undefined) => import("@mat3ra/esse/lib/js/types").EntityReferenceSchema);
        getEntityByName: ((entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string) => import("@exabyte-io/code.js/dist/entity").InMemoryEntity) & ((entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string) => import("@exabyte-io/code.js/dist/entity").InMemoryEntity) & ((entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string) => import("@exabyte-io/code.js/dist/entity").InMemoryEntity) & ((entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string) => import("@exabyte-io/code.js/dist/entity").InMemoryEntity) & ((entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string) => import("@exabyte-io/code.js/dist/entity").InMemoryEntity);
        metadata: object;
        updateMetadata(object: object): void;
        name: string;
        setName(name: string): void;
        readonly isDefault: boolean;
    };
    readonly jsonSchema: {
        $id: string;
        $schema: string;
        title: string;
        type: string;
        required: string[];
        properties: {
            formula: {
                description: string;
                type: string;
            };
            unitCellFormula: {
                description: string;
                type: string;
            };
            basis: {
                $schema: string;
                title: string;
                type: string;
                required: string[];
                properties: {
                    elements: {
                        type: string;
                        items: {
                            $schema: string;
                            title: string;
                            description: string;
                            type: string;
                            required: string[];
                            properties: {
                                id: {
                                    type: string;
                                };
                                value: {
                                    type: string;
                                };
                                occurrence: {
                                    description: string;
                                    type: string;
                                };
                                oxidationState: {
                                    type: string;
                                };
                            };
                        };
                    };
                    labels: {
                        description: string;
                        type: string;
                        items: {
                            properties: {
                                id: {
                                    type: string;
                                };
                                value: {
                                    type: string;
                                };
                            };
                        };
                    };
                    coordinates: {
                        type: string;
                        items: {
                            $schema: string;
                            title: string;
                            description: string;
                            type: string;
                            properties: {
                                id: {
                                    type: string;
                                };
                                value: {
                                    $schema: string;
                                    title: string;
                                    oneOf: {
                                        $schema: string;
                                        title: string;
                                        type: string;
                                        items: {
                                            type: string;
                                        };
                                        minItems: number;
                                        maxItems: number;
                                    }[];
                                };
                            };
                        };
                    };
                    name: {
                        type: string;
                    };
                    units: {
                        type: string;
                    };
                    bonds: {
                        $schema: string;
                        title: string;
                        type: string;
                        uniqueItems: boolean;
                        items: {
                            type: string;
                            properties: {
                                atomPair: {
                                    description: string;
                                    type: string;
                                    minItems: number;
                                    maxItems: number;
                                    $schema: string;
                                    title: string;
                                    items: {
                                        type: string;
                                        properties: {
                                            id: {
                                                description: string;
                                                type: string;
                                            };
                                        };
                                    };
                                };
                                bondType: {
                                    type: string;
                                    enum: string[];
                                };
                            }; /**
                             * Returns the inchi string from the derivedProperties for a non-periodic material, or throws an error if the
                             *  inchi cannot be found.
                             *  @returns {String}
                             */
                        };
                    };
                };
            };
            lattice: {
                $schema: string;
                title: string;
                type: string;
                required: string[];
                properties: {
                    name: {
                        enum: string[];
                    };
                    vectors: {
                        $schema: string;
                        title: string;
                        type: string;
                        required: string[];
                        properties: {
                            alat: {
                                description: string;
                                type: string;
                                default: number;
                            };
                            units: {
                                enum: string[];
                            };
                            a: {
                                $schema: string;
                                title: string;
                                type: string;
                                minItems: number;
                                maxItems: number;
                                items: {
                                    type: string;
                                };
                            };
                            b: {
                                $schema: string;
                                title: string;
                                type: string;
                                minItems: number;
                                maxItems: number;
                                items: {
                                    type: string;
                                };
                            };
                            c: {
                                $schema: string;
                                title: string;
                                type: string;
                                minItems: number;
                                maxItems: number;
                                items: {
                                    type: string;
                                };
                            };
                        };
                    };
                    type: {
                        $schema: string;
                        title: string;
                        type: string;
                        enum: string[];
                    };
                    units: {
                        type: string;
                        properties: {
                            length: {
                                type: string;
                                enum: string[];
                            };
                            angle: {
                                type: string;
                                enum: string[];
                            };
                        };
                    };
                    a: {
                        description: string;
                        type: string;
                    };
                    b: {
                        description: string;
                        type: string;
                    };
                    c: {
                        description: string;
                        type: string;
                    };
                    alpha: {
                        description: string;
                        type: string;
                    };
                    beta: {
                        description: string;
                        type: string;
                    };
                    gamma: {
                        description: string;
                        type: string;
                    };
                };
            };
            derivedProperties: {
                $schema: string;
                title: string;
                type: string;
                items: {
                    oneOf: ({
                        $schema: string;
                        title: string;
                        type: string;
                        required: string[];
                        properties: {
                            name: {
                                enum: string[];
                            };
                            units: {
                                enum: string[];
                            };
                            value: {
                                type: string;
                                minimum?: undefined;
                                maximum?: undefined;
                            };
                            pointGroupSymbol?: undefined;
                            spaceGroupSymbol?: undefined;
                            tolerance?: undefined;
                            element?: undefined;
                            degree?: undefined;
                        };
                        description?: undefined;
                    } | {
                        $schema: string;
                        title: string;
                        type: string;
                        properties: {
                            pointGroupSymbol: {
                                description: string;
                                type: string;
                            };
                            spaceGroupSymbol: {
                                description: string;
                                type: string;
                            };
                            tolerance: {
                                type: string;
                                description: string;
                                $schema: string;
                                title: string;
                                required: string[];
                                properties: {
                                    units: {
                                        enum: string[];
                                    };
                                    value: {
                                        type: string;
                                    };
                                };
                            };
                            name: {
                                enum: string[];
                            };
                            units?: undefined;
                            value?: undefined;
                            element?: undefined;
                            degree?: undefined;
                        };
                        required?: undefined;
                        description?: undefined;
                    } | {
                        $schema: string;
                        title: string;
                        description: string;
                        type: string;
                        required: string[];
                        properties: {
                            name: {
                                enum: string[];
                            };
                            value: {
                                type: string;
                                minimum: number;
                                maximum: number;
                            };
                            element: {
                                type: string;
                                description: string;
                            };
                            units?: undefined;
                            pointGroupSymbol?: undefined;
                            spaceGroupSymbol?: undefined;
                            tolerance?: undefined;
                            degree?: undefined;
                        };
                    } | {
                        $schema: string;
                        title: string;
                        description: string;
                        type: string;
                        required: string[];
                        properties: {
                            name: {
                                enum: string[];
                            };
                            degree: {
                                type: string;
                                description: string;
                            };
                            value: {
                                type: string;
                                minimum?: undefined;
                                maximum?: undefined;
                            };
                            units?: undefined;
                            pointGroupSymbol?: undefined;
                            spaceGroupSymbol?: undefined;
                            tolerance?: undefined;
                            element?: undefined;
                        };
                    } | {
                        $schema: string;
                        title: string;
                        type: string;
                        required: string[];
                        properties: {
                            name: {
                                enum: string[];
                            };
                            value: {
                                type: string;
                                minimum?: undefined;
                                maximum?: undefined;
                            };
                            units?: undefined;
                            pointGroupSymbol?: undefined;
                            spaceGroupSymbol?: undefined;
                            tolerance?: undefined;
                            element?: undefined;
                            degree?: undefined;
                        };
                        description?: undefined;
                    })[];
                    discriminator: {
                        propertyName: string;
                    };
                    required: string[];
                };
            };
            external: {
                $schema: string;
                title: string;
                description: string;
                type: string;
                required: string[];
                properties: {
                    id: {
                        description: string;
                        oneOf: {
                            type: string;
                        }[];
                    };
                    source: {
                        description: string;
                        type: string;
                    };
                    origin: {
                        description: string;
                        type: string;
                    };
                    data: {
                        description: string;
                        type: string;
                    };
                    doi: {
                        description: string;
                        type: string;
                    };
                    url: {
                        description: string;
                        type: string;
                    };
                };
            };
            src: {
                $schema: string;
                title: string;
                description: string;
                type: string;
                required: string[];
                properties: {
                    extension: {
                        description: string;
                        type: string;
                    };
                    filename: {
                        description: string;
                        type: string;
                    };
                    text: {
                        description: string;
                        type: string;
                    };
                    hash: {
                        description: string;
                        type: string;
                    };
                };
            };
            scaledHash: {
                description: string;
                type: string;
            };
            icsdId: {
                description: string;
                type: string;
            };
            isNonPeriodic: {
                description: string;
                type: string;
            };
            _id: {
                description: string;
                type: string;
            };
            slug: {
                description: string;
                type: string;
            };
            systemName: {
                type: string;
            };
            consistencyChecks: {
                type: string;
                items: {
                    $schema: string;
                    title: string;
                    type: string;
                    description: string;
                    required: string[];
                    properties: {
                        key: {
                            type: string;
                            description: string;
                        };
                        name: {
                            enum: string[];
                            description: string;
                        };
                        severity: {
                            enum: string[];
                            description: string;
                        };
                        message: {
                            type: string;
                            description: string;
                        };
                    };
                };
            };
            schemaVersion: {
                description: string;
                type: string;
                default: string;
            };
            name: {
                description: string;
                type: string;
            };
            isDefault: {
                description: string;
                type: string;
                default: boolean;
            };
            metadata: {
                type: string;
            };
        };
    };
    readonly defaultConfig: {
        name: string;
        basis: {
            elements: {
                id: number;
                value: string;
            }[];
            coordinates: {
                id: number;
                value: number[];
            }[];
            units: string;
        };
        lattice: {
            type: string;
            a: number;
            b: number;
            c: number;
            alpha: number;
            beta: number;
            gamma: number;
            units: {
                length: string;
                angle: string;
            };
        };
    };
} & (new (...args: any[]) => {
    consistencyChecks: object[];
    addConsistencyChecks(array: object[]): void;
    _json: AnyObject;
    prop<T = undefined>(name: string, defaultValue: T): T;
    prop<T_1 = undefined>(name: string): T_1 | undefined;
    setProp(name: string, value: unknown): void;
    unsetProp(name: string): void;
    toJSON(exclude?: string[] | undefined): AnyObject;
    toJSONSafe(exclude?: string[] | undefined): AnyObject;
    /**
     * Returns material's basis in XYZ format.
     */
    toJSONQuick(exclude?: string[] | undefined): AnyObject;
    clone(extraContext?: object | undefined): any;
    validate(): void;
    clean(config: AnyObject): AnyObject;
    isValid(): boolean;
    id: string;
    readonly cls: string;
    getClsName(): string;
    readonly slug: string;
    readonly isSystemEntity: boolean;
    getAsEntityReference(byIdOnly?: boolean | undefined): import("@mat3ra/esse/lib/js/types").EntityReferenceSchema;
    getEntityByName(entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string): import("@exabyte-io/code.js/dist/entity").InMemoryEntity;
}) & (new (...args: any[]) => {
    metadata: object;
    updateMetadata(object: object): void;
    _json: AnyObject;
    prop<T_2 = undefined>(name: string, defaultValue: T_2): T_2;
    prop<T_1_1 = undefined>(name: string): T_1_1 | undefined;
    setProp(name: string, value: unknown): void;
    unsetProp(name: string): void;
    toJSON(exclude?: string[] | undefined): AnyObject;
    toJSONSafe(exclude?: string[] | undefined): AnyObject;
    toJSONQuick(exclude?: string[] | undefined): AnyObject;
    clone(extraContext?: object | undefined): any;
    validate(): void;
    clean(config: AnyObject): AnyObject;
    isValid(): boolean;
    id: string;
    readonly cls: string;
    getClsName(): string;
    readonly slug: string;
    readonly isSystemEntity: boolean;
    getAsEntityReference(byIdOnly?: boolean | undefined): import("@mat3ra/esse/lib/js/types").EntityReferenceSchema;
    getEntityByName(entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string): import("@exabyte-io/code.js/dist/entity").InMemoryEntity;
}) & (new (...args: any[]) => {
    name: string;
    setName(name: string): void;
    _json: AnyObject;
    prop<T_3 = undefined>(name: string, defaultValue: T_3): T_3;
    prop<T_1_2 = undefined>(name: string): T_1_2 | undefined;
    setProp(name: string, value: unknown): void;
    unsetProp(name: string): void;
    toJSON(exclude?: string[] | undefined): AnyObject;
    toJSONSafe(exclude?: string[] | undefined): AnyObject;
    toJSONQuick(exclude?: string[] | undefined): AnyObject;
    clone(extraContext?: object | undefined): any;
    validate(): void;
    clean(config: AnyObject): AnyObject;
    isValid(): boolean;
    id: string;
    readonly cls: string;
    getClsName(): string;
    readonly slug: string;
    readonly isSystemEntity: boolean;
    /**
     * @summary a series of checks for the material's basis and returns an array of results in ConsistencyChecks format.
     * @returns Array of checks results
     */
    getAsEntityReference(byIdOnly?: boolean | undefined): import("@mat3ra/esse/lib/js/types").EntityReferenceSchema;
    getEntityByName(entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string): import("@exabyte-io/code.js/dist/entity").InMemoryEntity;
}) & {
    new (...args: any[]): {
        readonly isDefault: boolean;
        _json: AnyObject;
        prop<T_4 = undefined>(name: string, defaultValue: T_4): T_4;
        prop<T_1_3 = undefined>(name: string): T_1_3 | undefined;
        setProp(name: string, value: unknown): void;
        unsetProp(name: string): void;
        toJSON(exclude?: string[] | undefined): AnyObject;
        toJSONSafe(exclude?: string[] | undefined): AnyObject;
        toJSONQuick(exclude?: string[] | undefined): AnyObject;
        clone(extraContext?: object | undefined): any;
        validate(): void;
        clean(config: AnyObject): AnyObject;
        isValid(): boolean;
        id: string;
        readonly cls: string;
        getClsName(): string;
        readonly slug: string;
        readonly isSystemEntity: boolean;
        getAsEntityReference(byIdOnly?: boolean | undefined): import("@mat3ra/esse/lib/js/types").EntityReferenceSchema;
        getEntityByName(entities: import("@exabyte-io/code.js/dist/entity").InMemoryEntity[], entity: string, name: string): import("@exabyte-io/code.js/dist/entity").InMemoryEntity;
    };
    readonly defaultConfig: object | null;
    createDefault(): any;
} & typeof import("@exabyte-io/code.js/dist/entity").InMemoryEntity;
export type Material = InstanceType<typeof Material>;
export {};
