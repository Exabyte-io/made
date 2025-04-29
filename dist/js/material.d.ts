import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import { defaultMaterialConfig } from "./materialMixin";
export { defaultMaterialConfig };
declare const BaseEntity: (new (...args: any[]) => {
    consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
    addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
    _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    prop<T = undefined>(name: string, defaultValue: T): T;
    prop<T_1 = undefined>(name: string): T_1 | undefined;
    setProp(name: string, value: unknown): void;
    unsetProp(name: string): void;
    setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
    toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    clone(extraContext?: object | undefined): any;
    validate(): void;
    clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    isValid(): boolean;
    readonly cls: string;
    getClsName(): string;
    getAsEntityReference(byIdOnly: true): {
        _id: string;
    };
    getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
    getEntityByName(entities: import("@mat3ra/code/dist/js/entity").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity").InMemoryEntity;
    id: string;
    _id: string;
    schemaVersion: string;
    systemName: string;
    readonly slug: string;
    readonly isSystemEntity: boolean;
}) & (new (...args: any[]) => {
    metadata: object;
    updateMetadata(object: object): void;
    _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    prop<T_2 = undefined>(name: string, defaultValue: T_2): T_2;
    prop<T_3 = undefined>(name: string): T_3 | undefined;
    setProp(name: string, value: unknown): void;
    unsetProp(name: string): void;
    setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
    toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    clone(extraContext?: object | undefined): any;
    validate(): void;
    clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    isValid(): boolean;
    readonly cls: string;
    getClsName(): string;
    getAsEntityReference(byIdOnly: true): {
        _id: string;
    };
    getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
    getEntityByName(entities: import("@mat3ra/code/dist/js/entity").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity").InMemoryEntity;
    id: string;
    _id: string;
    schemaVersion: string;
    systemName: string;
    readonly slug: string;
    readonly isSystemEntity: boolean;
}) & typeof import("@mat3ra/code/dist/js/entity").InMemoryEntity & Constructor<{
    isDefault: boolean;
}> & {
    createDefault(this: Constructor<import("@mat3ra/code/dist/js/entity").InMemoryEntity> & {
        defaultConfig?: object | null | undefined;
    }): import("@mat3ra/code/dist/js/entity").InMemoryEntity;
} & import("@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin").NamedEntityConstructor;
type MaterialBaseEntity = Constructor<InstanceType<typeof BaseEntity>>;
export declare function MaterialMixin<T extends MaterialBaseEntity = MaterialBaseEntity>(superclass: T): {
    new (...config: any[]): {
        consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
        addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
        _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject & import("./materialMixin").MaterialSchemaJSON;
        prop: {
            <T_1 = undefined>(name: string, defaultValue: T_1): T_1;
            <T_2 = undefined>(name: string): T_2 | undefined;
        } & {
            <T_3 = undefined>(name: string, defaultValue: T_3): T_3;
            <T_4 = undefined>(name: string): T_4 | undefined;
        } & {
            <T_5 = undefined>(name: string, defaultValue: T_5): T_5;
            <T_6 = undefined>(name: string): T_6 | undefined;
        };
        setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
        unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
        setProps: ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any);
        toJSON: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & (() => import("./types").MaterialJSON);
        toJSONSafe: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
        toJSONQuick: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
        clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
        validate: (() => void) & (() => void) & (() => void);
        clean: ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
        isValid: (() => boolean) & (() => boolean) & (() => boolean);
        readonly cls: string;
        getClsName: (() => string) & (() => string) & (() => string);
        getAsEntityReference: {
            (byIdOnly: true): {
                _id: string;
            };
            (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
        } & {
            (byIdOnly: true): {
                _id: string;
            };
            (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
        } & {
            (byIdOnly: true): {
                _id: string;
            };
            (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
        };
        getEntityByName: ((entities: import("@mat3ra/code/dist/js/entity").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity").InMemoryEntity);
        id: string;
        _id: string;
        schemaVersion: string;
        systemName: string;
        readonly slug: string;
        readonly isSystemEntity: boolean;
        metadata: object;
        updateMetadata(object: object): void;
        isDefault: boolean;
        name: string;
        setName(name: string): void;
        src: import("@mat3ra/esse/dist/js/types").FileSourceSchema | undefined;
        isNonPeriodic: boolean;
        readonly formula: string;
        readonly unitCellFormula: string;
        readonly basis: import("./materialMixin").OptionallyConstrainedBasisConfig;
        readonly Basis: import("./basis/constrained_basis").ConstrainedBasis;
        readonly uniqueElements: string[];
        lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
        readonly Lattice: import("./lattice/lattice").Lattice;
        hash: string;
        readonly scaledHash: string;
        calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
        getInchiStringForHash(): string;
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
        getDerivedProperties(): import("@mat3ra/esse/dist/js/types").DerivedPropertiesSchema;
        unsetFileProps(): void;
        setBasis(textOrObject: string | import("./basis/basis").BasisConfig, format?: string | undefined, unitz?: string | undefined): void;
        updateFormula(): void;
        setBasisConstraints(constraints: import("./constraints/constraints").Constraint[]): void;
        toCrystal(): void;
        toCartesian(): void;
        getBasisAsXyz(fractional?: boolean): string;
        getAsQEFormat(): string;
        getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
        getACopyWithConventionalCell(): import("@mat3ra/code/dist/js/entity").InMemoryEntity & {
            name: string;
        } & {
            setName(name: string): void;
        };
        getConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
        getBasisConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
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
        constructMaterialFileSource(fileName: string, fileContent: string, fileExtension: string): import("@mat3ra/esse/dist/js/types").FileSourceSchema;
    };
} & T;
export declare const Material: {
    new (...config: any[]): {
        consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
        addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
        _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject & import("./materialMixin").MaterialSchemaJSON;
        prop: {
            <T = undefined>(name: string, defaultValue: T): T;
            <T_1 = undefined>(name: string): T_1 | undefined;
        } & {
            <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
            <T_3 = undefined>(name: string): T_3 | undefined;
        } & {
            <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
            <T_5 = undefined>(name: string): T_5 | undefined;
        };
        setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
        unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
        setProps: ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any);
        toJSON: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & (() => import("./types").MaterialJSON);
        toJSONSafe: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
        toJSONQuick: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
        clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
        validate: (() => void) & (() => void) & (() => void);
        clean: ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
        isValid: (() => boolean) & (() => boolean) & (() => boolean);
        readonly cls: string;
        getClsName: (() => string) & (() => string) & (() => string);
        getAsEntityReference: {
            (byIdOnly: true): {
                _id: string;
            };
            (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
        } & {
            (byIdOnly: true): {
                _id: string;
            };
            (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
        } & {
            (byIdOnly: true): {
                _id: string;
            };
            (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
        };
        getEntityByName: ((entities: import("@mat3ra/code/dist/js/entity").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity").InMemoryEntity);
        id: string;
        _id: string;
        schemaVersion: string;
        systemName: string;
        readonly slug: string;
        readonly isSystemEntity: boolean;
        metadata: object;
        updateMetadata(object: object): void;
        isDefault: boolean;
        name: string;
        setName(name: string): void;
        src: import("@mat3ra/esse/dist/js/types").FileSourceSchema | undefined;
        isNonPeriodic: boolean;
        readonly formula: string;
        readonly unitCellFormula: string;
        readonly basis: import("./materialMixin").OptionallyConstrainedBasisConfig;
        readonly Basis: import("./basis/constrained_basis").ConstrainedBasis;
        readonly uniqueElements: string[];
        lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
        readonly Lattice: import("./lattice/lattice").Lattice;
        hash: string;
        readonly scaledHash: string;
        calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
        getInchiStringForHash(): string;
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
        getDerivedProperties(): import("@mat3ra/esse/dist/js/types").DerivedPropertiesSchema;
        unsetFileProps(): void;
        setBasis(textOrObject: string | import("./basis/basis").BasisConfig, format?: string | undefined, unitz?: string | undefined): void;
        updateFormula(): void;
        setBasisConstraints(constraints: import("./constraints/constraints").Constraint[]): void;
        toCrystal(): void;
        toCartesian(): void;
        getBasisAsXyz(fractional?: boolean): string;
        getAsQEFormat(): string;
        getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
        getACopyWithConventionalCell(): import("@mat3ra/code/dist/js/entity").InMemoryEntity & {
            name: string;
        } & {
            setName(name: string): void;
        };
        getConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
        getBasisConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
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
        constructMaterialFileSource(fileName: string, fileContent: string, fileExtension: string): import("@mat3ra/esse/dist/js/types").FileSourceSchema;
    };
} & (new (...args: any[]) => {
    consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
    addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
    _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    prop<T = undefined>(name: string, defaultValue: T): T;
    prop<T_1 = undefined>(name: string): T_1 | undefined;
    setProp(name: string, value: unknown): void;
    unsetProp(name: string): void;
    setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
    toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    clone(extraContext?: object | undefined): any;
    validate(): void;
    clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    isValid(): boolean;
    readonly cls: string;
    getClsName(): string;
    getAsEntityReference(byIdOnly: true): {
        _id: string;
    };
    getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
    getEntityByName(entities: import("@mat3ra/code/dist/js/entity").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity").InMemoryEntity;
    id: string;
    _id: string;
    schemaVersion: string;
    systemName: string;
    readonly slug: string;
    readonly isSystemEntity: boolean;
}) & (new (...args: any[]) => {
    metadata: object;
    updateMetadata(object: object): void;
    _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    prop<T_2 = undefined>(name: string, defaultValue: T_2): T_2;
    prop<T_3 = undefined>(name: string): T_3 | undefined;
    setProp(name: string, value: unknown): void;
    unsetProp(name: string): void;
    setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
    toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    clone(extraContext?: object | undefined): any;
    validate(): void;
    clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
    isValid(): boolean;
    readonly cls: string;
    getClsName(): string;
    getAsEntityReference(byIdOnly: true): {
        _id: string;
    };
    getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
    getEntityByName(entities: import("@mat3ra/code/dist/js/entity").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity").InMemoryEntity;
    id: string;
    _id: string;
    schemaVersion: string;
    systemName: string;
    readonly slug: string;
    readonly isSystemEntity: boolean;
}) & typeof import("@mat3ra/code/dist/js/entity").InMemoryEntity & Constructor<{
    isDefault: boolean;
}> & {
    createDefault(this: Constructor<import("@mat3ra/code/dist/js/entity").InMemoryEntity> & {
        defaultConfig?: object | null | undefined;
    }): import("@mat3ra/code/dist/js/entity").InMemoryEntity;
} & import("@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin").NamedEntityConstructor;
export type Material = InstanceType<typeof Material>;
