"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Material = exports.defaultMaterialConfig = void 0;
const entity_1 = require("@mat3ra/code/dist/js/entity");
const materialMixin_1 = require("./materialMixin");
Object.defineProperty(exports, "defaultMaterialConfig", { enumerable: true, get: function () { return materialMixin_1.defaultMaterialConfig; } });
const BaseInMemoryEntity = entity_1.HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity;
class Material extends BaseInMemoryEntity {
    constructor(config) {
        super(config);
        (0, materialMixin_1.materialMixin)(this);
    }
}
exports.Material = Material;
(0, materialMixin_1.materialMixinStaticProps)(Material);
