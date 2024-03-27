import { ArrayOf3NumberElementsSchema, LatticeTypeSchema } from "@mat3ra/esse/dist/js/types";

export const DEFAULT_LATTICE_UNITS = {
    // by default lattice vectors shall be measured in angstrom, angles - in degrees
    length: {
        angstrom: "angstrom",
    },
    angle: {
        degree: "degree",
    },
};

export enum LatticeTypeExtended {
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
    TRI_1b = "TRI_1b",
}

export type Vector = ArrayOf3NumberElementsSchema;

export interface VectorsAsArray extends Array<Vector> {
    0: Vector;
    1: Vector;
    2: Vector;
}

interface LatticeTypeConfig {
    label: string;
    code: LatticeTypeSchema;
    editables: string[];
    editablesConventional: string[];
}

export const LATTICE_TYPE_CONFIGS: LatticeTypeConfig[] = [
    {
        label: "Simple Cubic",
        code: "CUB",
        // editables for primitive cell => WARNING: not tested
        editables: ["a"],
        // editables for conventional cell, taken from the publication above
        editablesConventional: ["a"],
    },
    {
        label: "Face-centered Cubic",
        code: "FCC",
        editables: ["a"],
        editablesConventional: ["a"],
    },
    {
        label: "Body-centered Cubic",
        code: "BCC",
        editables: ["a"],
        editablesConventional: ["a"],
    },
    {
        label: "Tetragonal",
        code: "TET",
        editables: ["a", "c"],
        editablesConventional: ["a", "c"],
    },
    {
        label: "Body-centered Tetragonal",
        code: "BCT",
        editables: ["a"],
        editablesConventional: ["a", "c"],
    },
    {
        label: "Orthorombic",
        code: "ORC",
        editables: ["a", "b", "c"],
        editablesConventional: ["a", "b", "c"],
    },
    {
        label: "Orthorombic Face-centered",
        code: "ORCF",
        editables: ["a", "b", "c"],
        editablesConventional: ["a", "b", "c"],
    },
    {
        label: "Orthorombic Body-centered",
        code: "ORCI",
        editables: ["a", "alpha", "gamma"],
        editablesConventional: ["a", "b", "c"],
    },
    {
        label: "Orthorombic Base-centered",
        code: "ORCC",
        editables: ["a", "c", "alpha"],
        editablesConventional: ["a", "b", "c"],
    },
    {
        label: "Hexagonal",
        code: "HEX",
        editables: ["a", "c"],
        editablesConventional: ["a", "c"],
    },
    {
        label: "Rhombohedral",
        code: "RHL",
        editables: ["a", "alpha"],
        editablesConventional: ["a", "alpha"],
    },
    {
        label: "Monoclinic",
        code: "MCL",
        editables: ["a", "b", "c", "alpha"],
        editablesConventional: ["a", "b", "c", "alpha"],
    },
    {
        label: "Monoclinic Base-centered",
        code: "MCLC",
        editables: ["a", "c", "alpha", "gamma"],
        editablesConventional: ["a", "b", "c", "alpha"],
    },
    {
        label: "Triclinic",
        code: "TRI",
        editables: ["a", "b", "c", "alpha", "beta", "gamma"],
        editablesConventional: ["a", "b", "c", "alpha", "beta", "gamma"],
    },
];
