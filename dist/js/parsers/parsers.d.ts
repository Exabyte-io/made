declare const _default: {
    xyz: {
        validate: typeof import("./xyz").validate;
        fromMaterial: (materialOrConfig: import("@mat3ra/esse/lib/js/types").MaterialSchema, fractional?: boolean) => string;
        toBasisConfig: (txt: string, units?: string, cell?: [import("@mat3ra/esse/lib/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/lib/js/types").ArrayOf3NumberElementsSchema, import("@mat3ra/esse/lib/js/types").ArrayOf3NumberElementsSchema]) => import("./xyz").BasisConfig;
        fromBasis: (basisClsInstance: import("../basis/constrained_basis").ConstrainedBasis, printFormat?: string, skipRounding?: boolean) => string;
        CombinatorialBasis: typeof import("./xyz_combinatorial_basis").CombinatorialBasis;
    };
    poscar: {
        isPoscar: (text: string) => boolean;
        toPoscar: (materialOrConfig: import("../types").MaterialJSON, omitConstraints?: boolean) => string;
        fromPoscar: (fileContent: string) => object;
        atomicConstraintsCharFromBool: (bool: boolean) => string;
        atomsCount: typeof import("./poscar").atomsCount;
    };
    cif: {
        parseMeta: (txt: string) => import("./cif").Meta;
    };
    espresso: {
        toEspressoFormat: (materialOrConfig: import("@mat3ra/esse/lib/js/types").MaterialSchema) => string;
    };
    nativeFormatParsers: {
        detectFormat: (text: string) => {
            JSON: string;
            POSCAR: string;
            CIF: string;
            PWX: string;
            XYZ: string;
            UNKNOWN: string;
        };
        convertFromNativeFormat: (text: string) => Object;
    };
};
export default _default;
