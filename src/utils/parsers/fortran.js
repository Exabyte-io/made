import { regex } from "./settings";

/**
 * Extracts pairs from a string data using provided regex pattern and type.
 * If isArray is set to true, treats the value as an array.
 *
 * @param {String} data - The string data to extract pairs from.
 * @param {RegExp} regexPattern - The regex pattern to use for extracting pairs.
 * @param {Function | NumberConstructor} type - The type of the value.
 * @param {Boolean} [isArray=false] - Whether to treat the value as an array.
 *
 * @returns {Array} The extracted pairs. Each pair is represented as an array,
 *                  where the first element is the key and the second element is the value.
 *                  If isArray is true, the value is an array where the first element
 *                  is the index of the Fortran array element and the second element is the value.
 * @throws {Error} If an invalid type is provided.
 */
function extractPairs(data, regexPattern, type, isArray) {
    return Array.from(data.matchAll(regexPattern)).map((match) => {
        switch (type) {
            case Number:
                return isArray
                    ? [match[1], parseInt(match[2], 10), parseFloat(match[3])]
                    : [match[1], parseFloat(match[2])];
            case String:
                return isArray
                    ? [match[1], parseInt(match[2], 10), match[3]]
                    : [match[1], match[2]];
            case Boolean:
                return isArray
                    ? [match[1], parseInt(match[2], 10), match[3] === "true"]
                    : [match[1], match[2] === "true"];
            default:
                throw new Error("Invalid type");
        }
    });
}

/**
 * @summary Extracts an array of the key value pairs from a Fortran namelist.
 * @param {String} data
 * @returns {Object[]}
 */
function extractKeyValuePairs(data) {
    const output = {};
    const numberPairs = extractPairs(data, regex.numberKeyValue, Number);
    const stringPairs = extractPairs(data, regex.stringKeyValue, String);
    const booleanPairs = extractPairs(data, regex.booleanKeyValue, Boolean);
    const numberArrayPairs = extractPairs(data, regex.numberArrayKeyValue, Number, true);
    const stringArrayPairs = extractPairs(data, regex.stringArrayKeyValue, String, true);
    const booleanArrayPairs = extractPairs(data, regex.booleanArrayKeyValue, Boolean, true);

    [...numberPairs, ...stringPairs, ...booleanPairs].forEach((pair) => {
        // eslint-disable-next-line prefer-destructuring
        output[pair[0]] = pair[1];
    });

    [numberArrayPairs, stringArrayPairs, booleanArrayPairs].forEach((arrayPairs) => {
        arrayPairs.forEach(([key, index, value]) => {
            if (!output[key]) output[key] = [];
            output[key][index - 1] = value; // use the index to put the value at the correct position
        });
    });

    return output;
}

/**
 * @summary Extracts namelist data from a string.
 * @param {String} text
 * @returns {Object}
 */
function extractNamelistData(text) {
    const namelistNameRegex = /^&(\w+)/gm;
    const matches = Array.from(text.matchAll(namelistNameRegex));
    const namelistNames = matches.map((match) => match[1].toLowerCase());
    const namelists = {};

    namelistNames.forEach((namelistName) => {
        const _regex = regex.namelists(namelistName);
        const data = text.match(_regex)[2];

        namelists[namelistName] = extractKeyValuePairs(data);
    });
    return namelists;
}

/** s
 * Parses Fortran namelists and cards data from a string.
 *
 * @summary Parses Fortran namelists and cards data from a QE input file string.
 * @param {String} text - The text to parse.
 * @throws {Error} If no namelist data is found in `text`.
 * @throws {Error} If no cards data is found in `text`.
 * @returns {Object} An object containing the parsed namelist and cards data. The exact structure of this object will depend on the structure of the namelist and cards data in `text`.
 */
export function parseFortranFile(text) {
    let output = {};
    try {
        output = extractNamelistData(text);
    } catch (err) {
        throw new Error("Incorrect fortran file");
    }

    if (!output) {
        throw new Error("No namelists found");
    }

    const match = text.match(regex.cards);
    if (!match || !match[1]) {
        throw new Error("No cards found");
    } else {
        // eslint-disable-next-line prefer-destructuring
        output.cards = match[1];
    }

    return output;
}
