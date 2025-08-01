import type { ConsistencyCheck, MaterialSchema } from "@mat3ra/esse/dist/js/types";
import { type MaterialMixinConstructor, defaultMaterialConfig } from "./materialMixin";
export { defaultMaterialConfig };
declare const BaseInMemoryEntity: typeof import("@mat3ra/code/dist/js/entity").InMemoryEntity & import("@mat3ra/code/dist/js/entity/mixins/DefaultableMixin").DefaultableInMemoryEntityConstructor & {
    createDefault<T extends import("@mat3ra/code/dist/js/utils/types").Constructor<import("@mat3ra/code/dist/js/entity").InMemoryEntity> & {
        defaultConfig?: object | null | undefined;
    }>(this: T): InstanceType<T> & {
        isDefault: boolean;
    };
} & import("@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin").NamedInMemoryEntityConstructor & import("@mat3ra/code/dist/js/entity/mixins/HasMetadataMixin").HasMetadataInMemoryEntityConstructor & import("@mat3ra/code/dist/js/entity/mixins/HasConsistencyChecksMixin").HasConsistencyChecksInMemoryEntityConstructor;
type BaseMaterial = MaterialMixinConstructor & typeof BaseInMemoryEntity;
type MaterialSchemaWithConsistencyChecksAsString = Omit<MaterialSchema, "consistencyChecks"> & {
    consistencyChecks?: ConsistencyCheck[];
};
declare const Material_base: BaseMaterial;
export declare class Material extends Material_base implements MaterialSchemaWithConsistencyChecksAsString {
    constructor(config: MaterialSchema);
}
