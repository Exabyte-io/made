import { HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

import {
    type MaterialMixinConstructor,
    defaultMaterialConfig,
    materialMixin,
} from "./materialMixin";

export { defaultMaterialConfig };

const BaseEntity = HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity;

type MaterialBaseEntity = Constructor<InstanceType<typeof BaseEntity>>;

export function MaterialMixin<T extends MaterialBaseEntity = MaterialBaseEntity>(superclass: T) {
    class MadeMaterial extends (superclass as T & MaterialMixinConstructor) {
        // TODO: add constraints (and other properties if needed) to ESSE MaterialSchema, then uncomment the line below to allow validation
        // During validation of the Material entity, properties absent in ESSE schema get deleted.
        // static readonly jsonSchema = MaterialJSONSchemaObject;

        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        constructor(...config: any[]) {
            super(...config);
            materialMixin(this);
            this.name = this.name || this.formula;
        }
    }

    return MadeMaterial;
}

export const Material = MaterialMixin(BaseEntity);

export type Material = InstanceType<typeof Material>;
