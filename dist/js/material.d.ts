import type { MaterialSchema } from "@mat3ra/esse/dist/js/types";
import { type MaterialMixinConstructor, defaultMaterialConfig } from "./materialMixin";
export { defaultMaterialConfig };
declare const BaseInMemoryEntity: (new (...args: any[]) => {
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
}) & typeof import("@mat3ra/code/dist/js/entity").InMemoryEntity & import("@mat3ra/code/dist/js/utils/types").Constructor<{
    isDefault: boolean;
}> & {
    createDefault<T_4 extends import("@mat3ra/code/dist/js/utils/types").Constructor<import("@mat3ra/code/dist/js/entity").InMemoryEntity> & {
        defaultConfig?: object | null | undefined;
    }>(this: T_4): InstanceType<T_4> & {
        isDefault: boolean;
    };
} & import("@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin").NamedEntityConstructor;
type BaseMaterial = typeof BaseInMemoryEntity & MaterialMixinConstructor;
declare const Material_base: BaseMaterial;
export declare class Material extends Material_base {
    constructor(config: MaterialSchema);
}
