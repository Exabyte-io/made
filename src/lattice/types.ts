export const DEFAULT_LATTICE_UNITS = {
    // by default lattice vectors shall be measured in angstrom, angles - in degrees
    length: {
        angstrom: "angstrom",
    },
    angle: {
        degree: "degree",
    },
};

/**
 * Shortlist of lattice types (according to [AFLOW](https://arxiv.org/abs/1004.2974))
 * Convention used is derived from:
 *   Setyawan, W., & Curtarolo, S. (2010). High-throughput electronic band structure calculations:
 *   Challenges and tools. Computational Materials Science, 49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010
 * Lattice parameters a, b, c are stored for CONVENTIONAL lattice, however the default unit cell is calculated for
 */

export enum LatticeType {
    CUB = "CUB",
    FCC = "FCC",
    BCC = "BCC",
    TET = "TET",
    BCT = "BCT",
    ORC = "ORC",
    ORCF = "ORCF",
    ORCI = "ORCI",
    ORCC = "ORCC",
    HEX = "HEX",
    RHL = "RHL",
    MCL = "MCL",
    MCLC = "MCLC",
    TRI = "TRI",
}

export const LATTICE_TYPE = {
    CUB: "CUB",
    FCC: "FCC",
    BCC: "BCC",
    TET: "TET",
    BCT: "BCT",
    ORC: "ORC",
    ORCF: "ORCF",
    ORCI: "ORCI",
    ORCC: "ORCC",
    HEX: "HEX",
    RHL: "RHL",
    MCL: "MCL",
    MCLC: "MCLC",
    TRI: "TRI",
};

export const LATTICE_TYPE_EXTENDED = {
    BCC: "BCC",
    BCT_1: "BCT-1",
    BCT_2: "BCT-2",
    CUB: "CUB",
    FCC: "FCC",
    HEX: "HEX",
    MCL: "MCL",
    MCLC_1: "MCLC-1",
    MCLC_2: "MCLC-2",
    MCLC_3: "MCLC-3",
    MCLC_4: "MCLC-4",
    MCLC_5: "MCLC-5",
    ORC: "ORC",
    ORCC: "ORCC",
    ORCF_1: "ORCF-1",
    ORCF_2: "ORCF-2",
    ORCF_3: "ORCF-3",
    ORCI: "ORCI",
    RHL_1: "RHL-1",
    RHL_2: "RHL-2",
    TET: "TET",
    TRI_1a: "TRI_1a",
    TRI_2a: "TRI_2a",
    TRI_1b: "TRI_1b",
};

export interface Vector extends Array<number> {
    0: number;
    1: number;
    2: number;
}

export interface VectorsAsArray extends Array<Vector> {
    0: Vector;
    1: Vector;
    2: Vector;
}

interface LatticeTypeConfig {
    label: string;
    code: LatticeType;
    editables: string[];
    editablesConventional: string[];
}

export const LATTICE_TYPE_CONFIGS: LatticeTypeConfig[] = [
    {
        label: "Simple Cubic",
        code: LatticeType.CUB,
        // editables for primitive cell => WARNING: not tested
        editables: ["a"],
        // editables for conventional cell, taken from the publication above
        editablesConventional: ["a"],
    },
    {
        label: "Face-centered Cubic",
        code: LatticeType.FCC,
        editables: ["a"],
        editablesConventional: ["a"],
    },
    {
        label: "Body-centered Cubic",
        code: LatticeType.BCC,
        editables: ["a"],
        editablesConventional: ["a"],
    },
    {
        label: "Tetragonal",
        code: LatticeType.TET,
        editables: ["a", "c"],
        editablesConventional: ["a", "c"],
    },
    {
        label: "Body-centered Tetragonal",
        code: LatticeType.BCT,
        editables: ["a"],
        editablesConventional: ["a", "c"],
    },
    {
        label: "Orthorombic",
        code: LatticeType.ORC,
        editables: ["a", "b", "c"],
        editablesConventional: ["a", "b", "c"],
    },
    {
        label: "Orthorombic Face-centered",
        code: LatticeType.ORCF,
        editables: ["a", "b", "c"],
        editablesConventional: ["a", "b", "c"],
    },
    {
        label: "Orthorombic Body-centered",
        code: LatticeType.ORCI,
        editables: ["a", "alpha", "gamma"],
        editablesConventional: ["a", "b", "c"],
    },
    {
        label: "Orthorombic Base-centered",
        code: LatticeType.ORCC,
        editables: ["a", "c", "alpha"],
        editablesConventional: ["a", "b", "c"],
    },
    {
        label: "Hexagonal",
        code: LatticeType.HEX,
        editables: ["a", "c"],
        editablesConventional: ["a", "c"],
    },
    {
        label: "Rhombohedral",
        code: LatticeType.RHL,
        editables: ["a", "alpha"],
        editablesConventional: ["a", "alpha"],
    },
    {
        label: "Monoclinic",
        code: LatticeType.MCL,
        editables: ["a", "b", "c", "alpha"],
        editablesConventional: ["a", "b", "c", "alpha"],
    },
    {
        label: "Monoclinic Base-centered",
        code: LatticeType.MCLC,
        editables: ["a", "c", "alpha", "gamma"],
        editablesConventional: ["a", "b", "c", "alpha"],
    },
    {
        label: "Triclinic",
        code: LatticeType.TRI,
        editables: ["a", "b", "c", "alpha", "beta", "gamma"],
        editablesConventional: ["a", "b", "c", "alpha", "beta", "gamma"],
    },
];
