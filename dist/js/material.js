"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Material = exports.defaultMaterialConfig = void 0;
const entity_1 = require("@mat3ra/code/dist/js/entity");
const DefaultableMixin_1 = require("@mat3ra/code/dist/js/entity/mixins/DefaultableMixin");
const HasConsistencyChecksMixin_1 = require("@mat3ra/code/dist/js/entity/mixins/HasConsistencyChecksMixin");
const HasMetadataMixin_1 = require("@mat3ra/code/dist/js/entity/mixins/HasMetadataMixin");
const NamedEntityMixin_1 = require("@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin");
const materialMixin_1 = require("./materialMixin");
Object.defineProperty(exports, "defaultMaterialConfig", { enumerable: true, get: function () { return materialMixin_1.defaultMaterialConfig; } });
class Material extends entity_1.InMemoryEntity {
}
exports.Material = Material;
(0, NamedEntityMixin_1.namedEntityMixin)(Material.prototype);
(0, DefaultableMixin_1.defaultableEntityMixin)(Material);
(0, HasConsistencyChecksMixin_1.hasConsistencyChecksMixin)(Material.prototype);
(0, HasMetadataMixin_1.hasMetadataMixin)(Material.prototype);
(0, materialMixin_1.materialMixin)(Material);
