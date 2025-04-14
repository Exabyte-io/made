declare namespace _default {
    export { generateConfig };
}
export default _default;
declare function generateConfig(material: any, millerIndices: any, numberOfLayers?: number, vx?: number, vy?: number): {
    name: string;
    basis: import("@mat3ra/esse/dist/js/types").BasisSchema;
    lattice: import("@mat3ra/esse/dist/js/types").LatticeImplicitSchema;
    outOfPlaneAxisIndex: any;
};
