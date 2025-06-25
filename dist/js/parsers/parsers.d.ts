declare const _default: {
    xyz: {
        validate: typeof import("./xyz").validate;
        fromMaterial: (materialOrConfig: import("@mat3ra/esse/dist/js/types").MaterialSchema, fractional?: boolean) => string;
        toBasisConfig: (txt: string, units?: string, cell?: import("../cell/cell").Cell) => import("../basis/constrained_basis").ConstrainedBasisConfig;
        fromBasis: (basisClsInstance: import("../basis/constrained_basis").ConstrainedBasis, coordinatePrintFormat: string) => string;
        CombinatorialBasis: typeof import("./xyz_combinatorial_basis").CombinatorialBasis;
    };
    poscar: {
        isPoscar: (text: string) => boolean;
        toPoscar: (materialOrConfig: import("../materialMixin").MaterialJSON, omitConstraints?: boolean) => string;
        fromPoscar: (fileContent: string) => object;
        atomicConstraintsCharFromBool: (bool: boolean) => string;
        atomsCount: typeof import("./poscar").atomsCount;
    };
    cif: {
        parseMeta: (txt: string) => import("./cif").Meta;
    };
    espresso: {
        toEspressoFormat: (materialOrConfig: import("@mat3ra/esse/dist/js/types").MaterialSchema) => string;
    };
    nativeFormatParsers: {
        detectFormat: (text: string) => string;
        convertFromNativeFormat: (text: string) => any;
    };
};
export default _default;
