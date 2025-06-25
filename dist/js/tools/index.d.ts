declare const _default: {
    surface: {
        generateConfig: (material: import("../material").Material, millerIndices: import("@mat3ra/esse/dist/js/types").Coordinate3DSchema, numberOfLayers?: number, vx?: number, vy?: number) => import("./surface").SlabConfigSchema;
    };
    supercell: {
        generateConfig: (material: import("../types").MaterialInMemoryEntity, supercellMatrix: import("@mat3ra/esse/dist/js/types").Matrix3X3Schema) => {
            name: string;
            basis: import("@mat3ra/esse/dist/js/types").BasisSchema;
            lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
        };
        generateNewBasisWithinSupercell: (basis: import("../types").Basis | import("../basis/constrained_basis").ConstrainedBasis, cell: import("../cell/cell").Cell, supercell: import("../cell/cell").Cell, supercellMatrix: import("@mat3ra/esse/dist/js/types").Matrix3X3Schema) => import("../types").Basis;
    };
    material: {
        scaleOneLatticeVector: (material: import("../material").Material, key?: "a" | "b" | "c", factor?: number) => void;
        scaleLatticeToMakeNonPeriodic: (material: import("../material").Material) => void;
        translateAtomsToCenter: (material: import("../material").Material) => void;
    };
    basis: {
        repeat: (basis: import("../types").Basis, repetitions: number[]) => import("../types").Basis;
        interpolate: (initialBasis: import("../types").Basis, finalBasis: import("../types").Basis, numberOfSteps?: number) => import("../types").Basis[];
    };
};
export default _default;
