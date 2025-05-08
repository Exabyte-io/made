declare const _default: {
    surface: {
        generateConfig: (material: import("../material").Material, millerIndices: import("@mat3ra/esse/dist/js/types").Coordinate3DSchema, numberOfLayers?: number, vx?: number, vy?: number) => import("./surface").SlabConfigSchema;
    };
    supercell: {
        generateConfig: (material: import("../types").MaterialInterface, supercellMatrix: import("@mat3ra/esse/dist/js/types").Matrix3X3Schema) => {
            name: string;
            basis: import("@mat3ra/esse/dist/js/types").BasisSchema;
            lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
        };
        generateNewBasisWithinSupercell: (basis: import("../basis/basis").Basis | import("../basis/constrained_basis").ConstrainedBasis, cell: import("../cell/cell").Cell, supercell: import("../cell/cell").Cell, supercellMatrix: import("@mat3ra/esse/dist/js/types").Matrix3X3Schema) => import("../basis/basis").Basis;
    };
    material: {
        scaleOneLatticeVector: (material: import("../material").Material, key?: "a" | "b" | "c", factor?: number) => void;
        scaleLatticeToMakeNonPeriodic: (material: import("../material").Material) => void;
        getBasisConfigTranslatedToCenter: (material: import("../material").Material) => void;
    };
    basis: {
        repeat: (basis: import("../basis/basis").Basis, repetitions: number[]) => import("../basis/basis").Basis;
        interpolate: (initialBasis: import("../basis/basis").Basis, finalBasis: import("../basis/basis").Basis, numberOfSteps?: number) => import("../basis/basis").Basis[];
    };
};
export default _default;
