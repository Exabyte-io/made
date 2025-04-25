import { Coordinate3DSchema, MaterialSchema } from "@mat3ra/esse/dist/js/types";
import { Material } from "../material";
export type SlabConfigSchema = MaterialSchema & {
    outOfPlaneAxisIndex: number;
};
declare function generateConfig(material: Material, millerIndices: Coordinate3DSchema, numberOfLayers?: number, vx?: number, vy?: number): SlabConfigSchema;
declare const _default: {
    generateConfig: typeof generateConfig;
};
export default _default;
