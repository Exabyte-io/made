import { HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { MaterialSchema } from "@mat3ra/esse/dist/js/types";

import {
    type MaterialMixinConstructor,
    defaultMaterialConfig,
    materialMixin,
    materialMixinStaticProps,
} from "./materialMixin";

export { defaultMaterialConfig };

const BaseInMemoryEntity = HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity;
type BaseMaterial = MaterialMixinConstructor & typeof BaseInMemoryEntity;

export class Material extends (BaseInMemoryEntity as BaseMaterial) {
    constructor(config: MaterialSchema) {
        super(config);
        materialMixin(this);
    }
}

materialMixinStaticProps(Material);
