import { LATTICE_TYPE } from "../../../lattice/types";
import { regex as fortranRegex } from "../../../utils/parsers/settings";

const double = fortranRegex.fortranDouble;
// Regexes for Espresso cards
export const regex = {
    espressoFingerprint: /&CONTROL|&SYSTEM|ATOMIC_SPECIES/i,
    atomicSpecies: new RegExp(
        "([A-Z][a-z]?)\\s" + // element Aa
            `${double}\\s` + // mass
            "(\\S*)\\s*" + // potential source file name
            "(?=\\n)", // end of line
        "gm",
    ),
    atomicPositionsUnits: new RegExp(
        "ATOMIC_POSITIONS\\s" + // start of card
            "\\(?" + // optional parentheses
            "(\\w+)" + // units
            "\\)?", // end of optional parentheses
    ),
    atomicPositions: new RegExp(
        `\\b([A-Z][a-z]*)\\b\\s+${double}\\s+${double}\\s+${double}(?:\\s+(0|1)\\s+(0|1)\\s+(0|1))?(?=\\s*\\n)`,
        "gm",
    ),
    cellParameters: new RegExp(
        `CELL_PARAMETERS\\s*(?:\\(?(?<units>\\w+)\\)?)?` +
            `\\s+(?<x1>${double})\\s+(?<y1>${double})\\s+(?<z1>${double})` +
            `\\s+(?<x2>${double})\\s+(?<y2>${double})\\s+(?<z2>${double})` +
            `\\s+(?<x3>${double})\\s+(?<y3>${double})\\s+(?<z3>${double})`,
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
