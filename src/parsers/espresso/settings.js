import { LATTICE_TYPE } from "../../lattice/types";
import { regex as commonRegex } from "../utils/settings";

const { double } = commonRegex.general;
export const regex = {
    espressoFingerprint: /&CONTROL|&SYSTEM|ATOMIC_SPECIES/i,
    atomicSpecies: new RegExp(
        "([A-Z][a-z]?)\\s+" + // element symbol Aa
            `(${double})\\s` + // mass
            "(\\S*)\\s*" + // potential source file name
            "(?=\\n)", // end of line
        "gm",
    ),
    atomicPositionsUnits: new RegExp(
        "ATOMIC_POSITIONS\\s+" + // start of card
            "\\(?" + // optional parentheses
            "(\\w+)" + // units
            "\\)?", // end of optional parentheses
    ),
    atomicPositions: new RegExp(
        `^\\s*([A-Z][a-z]*)\\s+` + // atomic element symbol
            `(${double})\\s+(${double})\\s+(${double})` + // atomic coordinates
            `(?:\\s+(0|1)\\s+(0|1)\\s+(0|1))?(?=\\s*\\n)`, // atomic constraints
        "gm",
    ),
    cellParameters: new RegExp(
        `CELL_PARAMETERS\\s*(?:\\(?(\\w+)\\)?)?\\n` +
            `^\\s*(${double})\\s+(${double})\\s+(${double})\\s*\\n` +
            `^\\s*(${double})\\s+(${double})\\s+(${double})\\s*\\n` +
            `^\\s*(${double})\\s+(${double})\\s+(${double})\\s*\\n`,
        "gm",
    ),
};

export const IBRAV_TO_LATTICE_TYPE_MAP = {
    1: LATTICE_TYPE.CUB,
    2: LATTICE_TYPE.FCC,
    3: LATTICE_TYPE.BCC,
    "-3": LATTICE_TYPE.BCC,
    4: LATTICE_TYPE.HEX,
    5: LATTICE_TYPE.RHL,
    "-5": LATTICE_TYPE.RHL,
    6: LATTICE_TYPE.TET,
    7: LATTICE_TYPE.BCT,
    8: LATTICE_TYPE.ORC,
    9: LATTICE_TYPE.ORCC,
    "-9": LATTICE_TYPE.ORCC,
    10: LATTICE_TYPE.ORCF,
    11: LATTICE_TYPE.ORCI,
    12: LATTICE_TYPE.MCL,
    "-12": LATTICE_TYPE.MCL,
    13: LATTICE_TYPE.MCLC,
    "-13": LATTICE_TYPE.MCLC,
    14: LATTICE_TYPE.TRI,
};
