import { regex as fortranRegex } from "../../utils/parsers/settings";

const QE_ENUMS = {
    controlNamelist: {
        header: "&CONTROL",
    },
    systemNamelist: {
        header: "&SYSTEM",
    },
    atomicSpeciesCard: {
        header: "ATOMIC_SPECIES",
    },
    atomicPositionsCard: {
        header: "ATOMIC_POSITIONS",
    },
    cellParametersCard: {
        header: "CELL_PARAMETERS",
    },
};

const double = fortranRegex.fortranDouble;
// Regexes for Espresso cards
export const regex = {
    espressoFingerprint: new RegExp(
        `${QE_ENUMS.controlNamelist.header}|${QE_ENUMS.systemNamelist.header}|${QE_ENUMS.atomicSpeciesCard.header}`,
        "i",
    ),
    atomicSpecies: new RegExp(`([A-Z][a-z]?)\\s+${double}\\s+(\\S*)\\s*(?=\\n)`, "gm"),
    atomicPositionsUnits: new RegExp(`${QE_ENUMS.atomicPositionsCard.header}\\s\\(?(\\w+)\\)?`),
    atomicPositions: new RegExp(
        `\\b([A-Z][a-z]*)\\b\\s+${double}\\s+${double}\\s+${double}(?:\\s+(0|1)\\s+(0|1)\\s+(0|1))?(?=\\s*\\n)`,
        "gm",
    ),
    cellParameters: new RegExp(
        `CELL_PARAMETERS\\s*(?:\\(?(\\w+)\\)?)?` +
            `\\s+${double}\\s+${double}\\s+${double}` +
            `\\s+${double}\\s+${double}\\s+${double}` +
            `\\s+${double}\\s+${double}\\s+${double}`,
    ),
};
