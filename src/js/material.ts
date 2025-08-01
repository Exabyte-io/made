import { HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { ConsistencyCheck, MaterialSchema } from "@mat3ra/esse/dist/js/types";

import {
    type MaterialMixinConstructor,
    defaultMaterialConfig,
    materialMixin,
    materialMixinStaticProps,
} from "./materialMixin";

export { defaultMaterialConfig };

const BaseInMemoryEntity = HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity;

type BaseMaterial = MaterialMixinConstructor & typeof BaseInMemoryEntity;

// TODO: remove in-line type creation
type MaterialSchemaWithConsistencyChecksAsString = Omit<MaterialSchema, "consistencyChecks"> & {
    consistencyChecks?: ConsistencyCheck[];
};

export class Material
    extends (BaseInMemoryEntity as BaseMaterial)
    implements MaterialSchemaWithConsistencyChecksAsString
{
    constructor(config: MaterialSchema) {
        super(config);
        materialMixin(this);
    }
}

materialMixinStaticProps(Material);
