import { regex as fortranRegex } from "../fortran/settings";
import { formatString } from "../utils";

// Regexes for Espresso cards
export const regex = {
    espressoFingerprint: /&CONTROL|&SYSTEM|ATOMIC_POSITIONS/i,
    atomicSpecies: new RegExp(
        formatString("([A-Z][a-z]?)\\s+{0}\\s+([\\S]*)\\s*(?=\\n)", fortranRegex.fortranDouble),
        "gm",
    ),
    atomicPositions: new RegExp(
        formatString(
            "\\b([A-Z][a-z]*)\\b\\s+{0}\\s+{0}\\s+{0}(?:\\s+(0|1)\\s+(0|1)\\s+(0|1))?(?=\\s*\\n)",
            fortranRegex.fortranDouble,
        ),
        "gm",
    ),
    atomicPositionsUnits: /ATOMIC_POSITIONS\s\(?(\w+)\)?/,
    cellParameters: new RegExp(
        formatString(
            "CELL_PARAMETERS\\s*(?:\\(?(\\w+)\\)?)?" +
                "\\s+{0}\\s+{0}\\s+{0}" +
                "\\s+{0}\\s+{0}\\s+{0}" +
                "\\s+{0}\\s+{0}\\s+{0}",
            fortranRegex.fortranDouble,
        ),
    ),
};
