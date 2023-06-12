const fortranDoubleRegex = "([-+]?\\d*\\.?\\d*(?:[Eed][+-]?\\d+)?)";
const fortranNamelistRegex = "&({0})((?:\\s|\\S)*?)\\/"; // capturing group: 1 - namelist name, 2 - data
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
function formatString(str, ...values) {
    return str.replace(/{(\d+)}/g, (match, index) => {
        return typeof values[index] !== "undefined" ? values[index] : match;
    });
}

// Regexes for the namelists & cards
const regex = {
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

/**
 * @summary Extracts an array of the key value pairs from a Fortran namelist.
 * @param {String} data
 * @returns {Object[]}
 */
function extractKeyValuePairs(data) {
    const output = {};

    const numberPairs = Array.from(data.matchAll(regex.numberKeyValue)).map((match) => [
        match[1],
        parseFloat(match[2]),
    ]);
    const stringPairs = Array.from(data.matchAll(regex.stringKeyValue)).map((match) => [
        match[1],
        match[2],
    ]);
    const booleanPairs = Array.from(data.matchAll(regex.booleanKeyValue)).map((match) => [
        match[1],
        match[2] === "true",
    ]);
    const numberArrayPairs = Array.from(data.matchAll(regex.numberArrayKeyValue)).map((match) => [
        match[1],
        parseFloat(match[3]),
    ]);

    [...numberPairs, ...stringPairs, ...booleanPairs].forEach((pair) => {
        // eslint-disable-next-line prefer-destructuring
        output[pair[0]] = pair[1];
    });

    numberArrayPairs.forEach((pair) => {
        if (!output[pair[0]]) output[pair[0]] = [null]; // to have celldm[i] start with 1 set the 0th element as null and push to it
        output[pair[0]].push(pair[1]);
    });

    return output;
}

/**
 * @summary Extracts namelist data from a string.
 * @param {String} text
 * @returns {Object}
 */
function extractNamelistData(text) {
    const namelistNameRegex = /&(\w+)/g;
    const matches = Array.from(text.matchAll(namelistNameRegex));
    const namelistNames = matches.map((match) => match[1].toLowerCase());

    // Create an object to hold all the key-value pairs for each namelist
    const namelists = {};

    // Iterate through each provided namelist name
    namelistNames.forEach((namelistName) => {
        // Create a new RegExp for the current namelist name
        const _regex = new RegExp(formatString(fortranNamelistRegex, namelistName.toUpperCase()));

        // Find the data for the current namelist
        const data = text.match(_regex)[2];

        // Extract the key-value pairs and store them in the namelists object
        namelists[namelistName] = extractKeyValuePairs(data);
    });
    return namelists;
}

/**
 * @summary Parses Fortran namelists and cards data from a string for a QE input file
 * @param {String} text
 * @returns {Object}
 */
export function parseFortranFile(text) {
    const output = extractNamelistData(text);

    // Cards are retrieved in a dedicated fashion:
    const cellParameters = text.match(regex.cellParameters);
    if (cellParameters) {
        output.cell = {
            cell: [
                [
                    parseFloat(cellParameters[2]),
                    parseFloat(cellParameters[3]),
                    parseFloat(cellParameters[4]),
                ],
                [
                    parseFloat(cellParameters[5]),
                    parseFloat(cellParameters[6]),
                    parseFloat(cellParameters[7]),
                ],
                [
                    parseFloat(cellParameters[8]),
                    parseFloat(cellParameters[9]),
                    parseFloat(cellParameters[10]),
                ],
            ],
            units: cellParameters[1],
        };
    }

    const atomicSpeciesMatches = Array.from(text.matchAll(regex.atomicSpecies));
    output.atomicSpecies = atomicSpeciesMatches.map((match) => ({
        element: match[1],
        mass: parseFloat(match[2]),
        potential: match[3],
    }));
    const atomicPositionsMatches = Array.from(text.matchAll(regex.atomicPositions));
    output.elements = atomicPositionsMatches.map((match, index) => ({
        id: index,
        value: match[1],
    }));
    output.coordinates = atomicPositionsMatches.map((match, index) => ({
        id: index,
        value: [parseFloat(match[2]), parseFloat(match[3]), parseFloat(match[4])],
    }));
    output.constraints = atomicPositionsMatches
        .filter((match) => match[5] && match[6] && match[7]) // Check if all three constrtaints exist
        .map((match, index) => ({
            id: index,
            value: [match[5] === "1", match[6] === "1", match[7] === "1"],
        }));
    // eslint-disable-next-line prefer-destructuring
    output.units = text.match(regex.atomicPositionsUnits)[1];

    return output;
}
