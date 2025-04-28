"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Material = exports.MaterialMixin = exports.defaultMaterialConfig = void 0;
const entity_1 = require("@mat3ra/code/dist/js/entity");
const materialMixin_1 = require("./materialMixin");
Object.defineProperty(exports, "defaultMaterialConfig", { enumerable: true, get: function () { return materialMixin_1.defaultMaterialConfig; } });
const BaseEntity = entity_1.HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity;
function MaterialMixin(superclass) {
    class MadeMaterial extends superclass {
        // TODO: add constraints (and other properties if needed) to ESSE MaterialSchema, then uncomment the line below to allow validation
        // During validation of the Material entity, properties absent in ESSE schema get deleted.
        // static readonly jsonSchema = MaterialJSONSchemaObject;
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        constructor(...config) {
            super(...config);
            (0, materialMixin_1.materialMixin)(this);
            this.name = this.name || this.formula;
        }
    }
    return MadeMaterial;
}
exports.MaterialMixin = MaterialMixin;
exports.Material = MaterialMixin(BaseEntity);
