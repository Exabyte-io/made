import { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import {
    type DefaultableInMemoryEntityConstructor,
    defaultableEntityMixin,
} from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import {
    HasConsistencyChecksInMemoryEntityConstructor,
    hasConsistencyChecksMixin,
} from "@mat3ra/code/dist/js/entity/mixins/HasConsistencyChecksMixin";
import {
    type HasMetadataInMemoryEntityConstructor,
    hasMetadataMixin,
} from "@mat3ra/code/dist/js/entity/mixins/HasMetadataMixin";
import {
    type NamedInMemoryEntityConstructor,
    namedEntityMixin,
} from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { ConsistencyCheck, MaterialSchema } from "@mat3ra/esse/dist/js/types";

import {
    type MaterialMixinConstructor,
    defaultMaterialConfig,
    materialMixin,
} from "./materialMixin";

export { defaultMaterialConfig };

type BaseMaterial = typeof InMemoryEntity &
    HasConsistencyChecksInMemoryEntityConstructor &
    DefaultableInMemoryEntityConstructor &
    HasMetadataInMemoryEntityConstructor &
    NamedInMemoryEntityConstructor &
    MaterialMixinConstructor;

// TODO: remove in-line type creation
type MaterialSchemaWithConsistencyChecksAsString = Omit<MaterialSchema, "consistencyChecks"> & {
    consistencyChecks?: ConsistencyCheck[];
};

type Schema = MaterialSchemaWithConsistencyChecksAsString;

export class Material extends (InMemoryEntity as BaseMaterial) implements Schema {}

namedEntityMixin(Material.prototype);
defaultableEntityMixin(Material);
hasConsistencyChecksMixin(Material.prototype);
hasMetadataMixin(Material.prototype);
materialMixin(Material);
