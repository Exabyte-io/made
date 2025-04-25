declare const _default: {
    surface: {
        generateConfig: (material: {
            _json: import("../material").MaterialSchemaJSON;
            toJSON(): import("../types/material").MaterialJSON;
            src: import("@mat3ra/esse/dist/js/types").FileSourceSchema | undefined;
            updateFormula(): void;
            isNonPeriodic: boolean;
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
            readonly formula: string;
            readonly unitCellFormula: string;
            unsetFileProps(): void;
            setBasis(textOrObject: string | import("../basis/basis").BasisConfig, format?: string | undefined, unitz?: string | undefined): void;
            setBasisConstraints(constraints: import("../constraints/constraints").Constraint[]): void;
            readonly basis: import("../material").OptionallyConstrainedBasisConfig;
            readonly Basis: import("../basis/constrained_basis").ConstrainedBasis;
            readonly uniqueElements: string[];
            lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
            readonly Lattice: import("../lattice/lattice").Lattice;
            getInchiStringForHash(): string;
            calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
            hash: string;
            readonly scaledHash: string;
            toCrystal(): void;
            toCartesian(): void;
            getBasisAsXyz(fractional?: boolean): string;
            getAsQEFormat(): string;
            getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
            getACopyWithConventionalCell(): any;
            getConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            getBasisConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
            prop: {
                <T = undefined>(name: string, defaultValue: T): T;
                <T_1 = undefined>(name: string): T_1 | undefined;
            } & {
                <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                <T_3 = undefined>(name: string): T_3 | undefined;
            } & {
                <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
                <T_5 = undefined>(name: string): T_5 | undefined;
            } & {
                <T_6 = undefined>(name: string, defaultValue: T_6): T_6;
                <T_7 = undefined>(name: string): T_7 | undefined;
            } & {
                <T_8 = undefined>(name: string, defaultValue: T_8): T_8;
                <T_9 = undefined>(name: string): T_9 | undefined;
            };
            setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
            unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
            setProps: ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any);
            toJSONSafe: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            toJSONQuick: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
            validate: (() => void) & (() => void) & (() => void) & (() => void) & (() => void);
            clean: ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            isValid: (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean);
            readonly cls: string;
            getClsName: (() => string) & (() => string) & (() => string) & (() => string) & (() => string);
            getAsEntityReference: {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            };
            getEntityByName: ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity);
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
            metadata: object;
            updateMetadata(object: object): void;
            name: string;
            setName(name: string): void;
            readonly isDefault: boolean;
        } & {
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
            name: string;
            setName(name: string): void;
            _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
            prop<T_4 = undefined>(name: string, defaultValue: T_4): T_4;
            prop<T_5 = undefined>(name: string): T_5 | undefined;
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
            readonly isDefault: boolean;
            _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
            prop<T_6 = undefined>(name: string, defaultValue: T_6): T_6;
            prop<T_7 = undefined>(name: string): T_7 | undefined;
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity, millerIndices: import("@mat3ra/esse/dist/js/types").Coordinate3DSchema, numberOfLayers?: number, vx?: number, vy?: number) => import("./surface").SlabConfigSchema;
    };
    supercell: {
        generateConfig: (material: import("../types").MaterialInterface, supercellMatrix: import("@mat3ra/esse/dist/js/types").Matrix3X3Schema) => {
            name: string;
            basis: import("@mat3ra/esse/dist/js/types").BasisSchema;
            lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
        };
        generateNewBasisWithinSupercell: (basis: import("../basis/basis").Basis | import("../basis/constrained_basis").ConstrainedBasis, cell: import("../cell/cell").Cell, supercell: import("../cell/cell").Cell, supercellMatrix: import("@mat3ra/esse/dist/js/types").Matrix3X3Schema) => import("../basis/basis").Basis;
    };
    material: {
        scaleOneLatticeVector: (material: {
            _json: import("../material").MaterialSchemaJSON;
            toJSON(): import("../types/material").MaterialJSON;
            src: import("@mat3ra/esse/dist/js/types").FileSourceSchema | undefined;
            updateFormula(): void;
            isNonPeriodic: boolean;
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
            readonly formula: string;
            readonly unitCellFormula: string;
            unsetFileProps(): void;
            setBasis(textOrObject: string | import("../basis/basis").BasisConfig, format?: string | undefined, unitz?: string | undefined): void;
            setBasisConstraints(constraints: import("../constraints/constraints").Constraint[]): void;
            readonly basis: import("../material").OptionallyConstrainedBasisConfig;
            readonly Basis: import("../basis/constrained_basis").ConstrainedBasis;
            readonly uniqueElements: string[];
            lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
            readonly Lattice: import("../lattice/lattice").Lattice;
            getInchiStringForHash(): string;
            calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
            hash: string;
            readonly scaledHash: string;
            toCrystal(): void;
            toCartesian(): void;
            getBasisAsXyz(fractional?: boolean): string;
            getAsQEFormat(): string;
            getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
            getACopyWithConventionalCell(): any;
            getConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            getBasisConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
            prop: {
                <T = undefined>(name: string, defaultValue: T): T;
                <T_1 = undefined>(name: string): T_1 | undefined;
            } & {
                <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                <T_3 = undefined>(name: string): T_3 | undefined;
            } & {
                <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
                <T_5 = undefined>(name: string): T_5 | undefined;
            } & {
                <T_6 = undefined>(name: string, defaultValue: T_6): T_6;
                <T_7 = undefined>(name: string): T_7 | undefined;
            } & {
                <T_8 = undefined>(name: string, defaultValue: T_8): T_8;
                <T_9 = undefined>(name: string): T_9 | undefined;
            };
            setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
            unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
            setProps: ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any);
            toJSONSafe: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            toJSONQuick: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
            validate: (() => void) & (() => void) & (() => void) & (() => void) & (() => void);
            clean: ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            isValid: (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean);
            readonly cls: string;
            getClsName: (() => string) & (() => string) & (() => string) & (() => string) & (() => string);
            getAsEntityReference: {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            };
            getEntityByName: ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity);
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
            metadata: object;
            updateMetadata(object: object): void;
            name: string;
            setName(name: string): void;
            readonly isDefault: boolean;
        } & {
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
            name: string;
            setName(name: string): void;
            _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
            prop<T_4 = undefined>(name: string, defaultValue: T_4): T_4;
            prop<T_5 = undefined>(name: string): T_5 | undefined;
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
            readonly isDefault: boolean;
            _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
            prop<T_6 = undefined>(name: string, defaultValue: T_6): T_6;
            prop<T_7 = undefined>(name: string): T_7 | undefined;
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity, key?: "a" | "b" | "c", factor?: number) => void;
        scaleLatticeToMakeNonPeriodic: (material: {
            _json: import("../material").MaterialSchemaJSON;
            toJSON(): import("../types/material").MaterialJSON;
            src: import("@mat3ra/esse/dist/js/types").FileSourceSchema | undefined;
            updateFormula(): void;
            isNonPeriodic: boolean;
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
            readonly formula: string;
            readonly unitCellFormula: string;
            unsetFileProps(): void;
            setBasis(textOrObject: string | import("../basis/basis").BasisConfig, format?: string | undefined, unitz?: string | undefined): void;
            setBasisConstraints(constraints: import("../constraints/constraints").Constraint[]): void;
            readonly basis: import("../material").OptionallyConstrainedBasisConfig;
            readonly Basis: import("../basis/constrained_basis").ConstrainedBasis;
            readonly uniqueElements: string[];
            lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
            readonly Lattice: import("../lattice/lattice").Lattice;
            getInchiStringForHash(): string;
            calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
            hash: string;
            readonly scaledHash: string;
            toCrystal(): void;
            toCartesian(): void;
            getBasisAsXyz(fractional?: boolean): string;
            getAsQEFormat(): string;
            getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
            getACopyWithConventionalCell(): any;
            getConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            getBasisConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
            prop: {
                <T = undefined>(name: string, defaultValue: T): T;
                <T_1 = undefined>(name: string): T_1 | undefined;
            } & {
                <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                <T_3 = undefined>(name: string): T_3 | undefined;
            } & {
                <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
                <T_5 = undefined>(name: string): T_5 | undefined;
            } & {
                <T_6 = undefined>(name: string, defaultValue: T_6): T_6;
                <T_7 = undefined>(name: string): T_7 | undefined;
            } & {
                <T_8 = undefined>(name: string, defaultValue: T_8): T_8;
                <T_9 = undefined>(name: string): T_9 | undefined;
            };
            setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
            unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
            setProps: ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any);
            toJSONSafe: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            toJSONQuick: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
            validate: (() => void) & (() => void) & (() => void) & (() => void) & (() => void);
            clean: ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            isValid: (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean);
            readonly cls: string;
            getClsName: (() => string) & (() => string) & (() => string) & (() => string) & (() => string);
            getAsEntityReference: {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            };
            getEntityByName: ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity);
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
            metadata: object;
            updateMetadata(object: object): void;
            name: string;
            setName(name: string): void;
            readonly isDefault: boolean;
        } & {
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
            name: string;
            setName(name: string): void;
            _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
            prop<T_4 = undefined>(name: string, defaultValue: T_4): T_4;
            prop<T_5 = undefined>(name: string): T_5 | undefined;
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
            readonly isDefault: boolean;
            _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
            prop<T_6 = undefined>(name: string, defaultValue: T_6): T_6;
            prop<T_7 = undefined>(name: string): T_7 | undefined;
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) => void;
        getBasisConfigTranslatedToCenter: (material: {
            _json: import("../material").MaterialSchemaJSON;
            toJSON(): import("../types/material").MaterialJSON;
            src: import("@mat3ra/esse/dist/js/types").FileSourceSchema | undefined;
            updateFormula(): void;
            isNonPeriodic: boolean;
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
            readonly formula: string;
            readonly unitCellFormula: string;
            unsetFileProps(): void;
            setBasis(textOrObject: string | import("../basis/basis").BasisConfig, format?: string | undefined, unitz?: string | undefined): void;
            setBasisConstraints(constraints: import("../constraints/constraints").Constraint[]): void;
            readonly basis: import("../material").OptionallyConstrainedBasisConfig;
            readonly Basis: import("../basis/constrained_basis").ConstrainedBasis;
            readonly uniqueElements: string[];
            lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
            readonly Lattice: import("../lattice/lattice").Lattice;
            getInchiStringForHash(): string;
            calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
            hash: string;
            readonly scaledHash: string;
            toCrystal(): void;
            toCartesian(): void;
            getBasisAsXyz(fractional?: boolean): string;
            getAsQEFormat(): string;
            getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
            getACopyWithConventionalCell(): any;
            getConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            getBasisConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
            prop: {
                <T = undefined>(name: string, defaultValue: T): T;
                <T_1 = undefined>(name: string): T_1 | undefined;
            } & {
                <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                <T_3 = undefined>(name: string): T_3 | undefined;
            } & {
                <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
                <T_5 = undefined>(name: string): T_5 | undefined;
            } & {
                <T_6 = undefined>(name: string, defaultValue: T_6): T_6;
                <T_7 = undefined>(name: string): T_7 | undefined;
            } & {
                <T_8 = undefined>(name: string, defaultValue: T_8): T_8;
                <T_9 = undefined>(name: string): T_9 | undefined;
            };
            setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
            unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
            setProps: ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any);
            toJSONSafe: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            toJSONQuick: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
            validate: (() => void) & (() => void) & (() => void) & (() => void) & (() => void);
            clean: ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            isValid: (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean) & (() => boolean);
            readonly cls: string;
            getClsName: (() => string) & (() => string) & (() => string) & (() => string) & (() => string);
            getAsEntityReference: {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            };
            getEntityByName: ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity);
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
            metadata: object;
            updateMetadata(object: object): void;
            name: string;
            setName(name: string): void;
            readonly isDefault: boolean;
        } & {
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
            name: string;
            setName(name: string): void;
            _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
            prop<T_4 = undefined>(name: string, defaultValue: T_4): T_4;
            prop<T_5 = undefined>(name: string): T_5 | undefined;
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & {
            readonly isDefault: boolean;
            _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
            prop<T_6 = undefined>(name: string, defaultValue: T_6): T_6;
            prop<T_7 = undefined>(name: string): T_7 | undefined;
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
            getAsEntityReference(byIdOnly: false): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
        } & import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) => void;
    };
    basis: {
        repeat: (basis: import("../basis/basis").Basis, repetitions: number[]) => import("../basis/basis").Basis;
        interpolate: (initialBasis: import("../basis/basis").Basis, finalBasis: import("../basis/basis").Basis, numberOfSteps?: number) => import("../basis/basis").Basis[];
    };
};
export default _default;
