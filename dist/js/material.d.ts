import { HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { ConsistencyCheck, MaterialSchema } from "@mat3ra/esse/dist/js/types";
import { type MaterialMixinConstructor, defaultMaterialConfig } from "./materialMixin";
export { defaultMaterialConfig };
declare const BaseInMemoryEntity: typeof HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity;
type BaseMaterial = MaterialMixinConstructor & typeof BaseInMemoryEntity;
type MaterialSchemaWithConsistencyChecksAsString = Omit<MaterialSchema, "consistencyChecks"> & {
    consistencyChecks?: ConsistencyCheck[];
};
type Schema = MaterialSchemaWithConsistencyChecksAsString;
declare const Material_base: BaseMaterial;
export declare class Material extends Material_base implements Schema {
    constructor(config: MaterialSchema);
}
