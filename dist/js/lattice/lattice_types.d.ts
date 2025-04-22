import { LatticeTypeEnum } from "@mat3ra/esse/dist/js/types";
export declare enum LatticeTypeExtended {
    BCC = "BCC",
    BCT_1 = "BCT-1",
    BCT_2 = "BCT-2",
    CUB = "CUB",
    FCC = "FCC",
    HEX = "HEX",
    MCL = "MCL",
    MCLC_1 = "MCLC-1",
    MCLC_2 = "MCLC-2",
    MCLC_3 = "MCLC-3",
    MCLC_4 = "MCLC-4",
    MCLC_5 = "MCLC-5",
    ORC = "ORC",
    ORCC = "ORCC",
    ORCF_1 = "ORCF-1",
    ORCF_2 = "ORCF-2",
    ORCF_3 = "ORCF-3",
    ORCI = "ORCI",
    RHL_1 = "RHL-1",
    RHL_2 = "RHL-2",
    TET = "TET",
    TRI_1a = "TRI_1a",
    TRI_2a = "TRI_2a",
    TRI_1b = "TRI_1b"
}
export interface LatticeTypeConfig {
    label: string;
    code: LatticeTypeEnum;
    editables: string[];
    editablesConventional: string[];
}
export declare const DEFAULT_LATTICE_UNITS: {
    length: {
        angstrom: string;
    };
    angle: {
        degree: string;
    };
};
export declare const LATTICE_TYPE_CONFIGS: LatticeTypeConfig[];
