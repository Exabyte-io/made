import { regex } from "./settings";

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
        parseInt(match[2], 10), // get the index of Fortran array element
        parseFloat(match[3]), // get the value of the Fortran array element
    ]);

    [...numberPairs, ...stringPairs, ...booleanPairs].forEach((pair) => {
        // eslint-disable-next-line prefer-destructuring
        output[pair[0]] = pair[1];
    });

    numberArrayPairs.forEach((pair) => {
        if (!output[pair[0]]) output[pair[0]] = [];
        // eslint-disable-next-line prefer-destructuring
        output[pair[0]][pair[1] - 1] = pair[2]; // use the index to put the value at the correct position
    });
    // TODO: add string and boolean arrays
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
