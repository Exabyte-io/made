const fortranDoubleRegex = "([-+]?\\d*\\.?\\d*(?:[Eed][+-]?\\d+)?)";
const fortranNamelistRegex = "&({0})((?:\\s|\\S)*?)\\/"; // capturing group: 1 - namelist name, 2 - data
const fortranCardsRegex = "^\\s*\\/(?![\\s\\S]*^\\/)([\\s\\S]*)"; // capturing group: 1 cards section string
const keyValueRegex = "^\\s*{0}\\s*=\\s*{1}\\s*\\n"; // 0 = string value name, 1 = str|bool|number
const stringRegex = "([+\\w.\\-\\/]*)"; // capturing group: 1 - string value
const fortranStringRegex = "'([\\w.\\-\\+\\/ ]*)'"; // capturing group: 1 - string value
const fortranArrayRegex = "^\\s*{0}\\({1}\\)\\s*=\\s*{2}\\s*\\n"; // /gm, capturing groups: 0 - array name, 1 - array index, 2 - value
const fortranBooleanRegex = "\\.(true|false)\\."; // capturing group: 1 - boolean value

/**
 * @summary Formats a string with other string to be used as a regex.
 * @param {String} str
 * @param {String | String[]} values
 * @returns {String}
 */
export function formatString(str, ...values) {
    return str.replace(/{(\d+)}/g, (match, index) => {
        return typeof values[index] !== "undefined" ? values[index] : match;
    });
}

// Regexes for the namelists & cards
export const regex = {
    stringKeyValue: new RegExp(formatString(keyValueRegex, stringRegex, fortranStringRegex), "gm"),
    numberKeyValue: new RegExp(formatString(keyValueRegex, stringRegex, fortranDoubleRegex), "gm"),
    booleanKeyValue: new RegExp(
        formatString(keyValueRegex, stringRegex, fortranBooleanRegex),
        "gm",
    ),
    numberArrayKeyValue: new RegExp(
        formatString(fortranArrayRegex, stringRegex, "(\\d+)", fortranDoubleRegex),
        "gm",
    ),
    namelists: (namelistName) =>
        new RegExp(formatString(fortranNamelistRegex, namelistName.toUpperCase())),
    cards: new RegExp(fortranCardsRegex, "m"),
    atomicSpecies: new RegExp(
        formatString("([A-Z][a-z]?)\\s+{0}\\s+([\\S]*)\\s*(?=\\n)", fortranDoubleRegex),
        "gm",
    ),
    atomicPositions: new RegExp(
        formatString(
            "\\b([A-Z][a-z]*)\\b\\s+{0}\\s+{0}\\s+{0}(?:\\s+(0|1)\\s+(0|1)\\s+(0|1))?(?=\\s*\\n)",
            fortranDoubleRegex,
        ),
        "gm",
    ),
    atomicPositionsUnits: /ATOMIC_POSITIONS\s\(?([\w]+)\)?/,
    cellParameters: new RegExp(
        formatString(
            "CELL_PARAMETERS\\s*(?:\\((\\w+)\\))?" +
                "\\s+{0}\\s+{0}\\s+{0}" +
                "\\s+{0}\\s+{0}\\s+{0}" +
                "\\s+{0}\\s+{0}\\s+{0}",
            fortranDoubleRegex,
        ),
    ),
};
