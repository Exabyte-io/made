import { HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { MaterialSchema } from "@mat3ra/esse/dist/js/types";
import { type MaterialMixinConstructor, defaultMaterialConfig } from "./materialMixin";
export { defaultMaterialConfig };
declare const BaseInMemoryEntity: typeof HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity;
type BaseMaterial = MaterialMixinConstructor & typeof BaseInMemoryEntity;
declare const Material_base: BaseMaterial;
export declare class Material extends Material_base {
    constructor(config: MaterialSchema);
}
