import { regex } from "./settings";

const typeParsers = {
    [Number]: (value) => parseFloat(value),
    [String]: (value) => value,
    [Boolean]: (value) => value === "true",
};

export const FortranParserMixin = (superclass) =>
    class extends superclass {
        /**
         * Parses Fortran namelists and cards data from a string.
         *
         * @summary Parses Fortran namelists and cards data from a QE input file string.
         * @param {String} content - The text to parse.
         * @throws {Error} If no namelist data is found in `text`.
         * @throws {Error} If no cards data is found in `text`.
         * @returns {Object} An object containing the parsed namelist and cards data. The exact structure of this object will depend on the structure of the namelist and cards data in `text`.
         */
        fortranParseNamelists(content) {
            let output = {};
            try {
                output = this.fortranExtractNamelistData(content);
            } catch (err) {
                throw new Error("Incorrect fortran file");
            }

            const match = regex.fortran.cards.exec(content);
            // eslint-disable-next-line prefer-destructuring
            output.cards = match[0];
            return output;
        }

        /**
         * @summary Extracts namelist data from a string.
         * @param {String} text
         * @returns {Object}
         */
        fortranExtractNamelistData(text) {
            const namelistNameRegex = /^&(\w+)/gm;
            const matches = Array.from(text.matchAll(namelistNameRegex));
            const namelistNames = matches.map((match) => match[1].toLowerCase());
            const namelists = {};

            namelistNames.forEach((namelistName) => {
                const namelistsRegex = regex.fortran.namelists(namelistName);
                const namelistData = text.match(namelistsRegex)[2];
                namelists[namelistName] = this.fortranGetKeyValuePairs(namelistData);
            });
            return namelists;
        }

        /**
         * @summary Extracts an array of the key value pairs from a Fortran namelist.
         * @param {String} data - namelist data
         * @returns {Object[]}
         */
        fortranGetKeyValuePairs(data) {
            const output = {};
            // TODO: Fix to convert numbers like 1.234D-567
            const numberPairs = this.fortranExtractKeyValuePairs(
                data,
                regex.fortran.numberKeyValue,
                Number,
            );
            const stringPairs = this.fortranExtractKeyValuePairs(
                data,
                regex.fortran.stringKeyValue,
                String,
            );
            const booleanPairs = this.fortranExtractKeyValuePairs(
                data,
                regex.fortran.booleanKeyValue,
                Boolean,
            );
            // TODO: Add functionality to parse Fortran lists assigned with multiple values inline: list = 1,2,3 -- current implementation doesn't capture that
            const numberArrayPairs = this.fortranExtractKeyValuePairs(
                data,
                regex.fortran.numberArrayKeyValue,
                Number,
                true,
            );
            const stringArrayPairs = this.fortranExtractKeyValuePairs(
                data,
                regex.fortran.stringArrayKeyValue,
                String,
                true,
            );
            const booleanArrayPairs = this.fortranExtractKeyValuePairs(
                data,
                regex.fortran.booleanArrayKeyValue,
                Boolean,
                true,
            );

            [...numberPairs, ...stringPairs, ...booleanPairs].forEach((pair) => {
                // eslint-disable-next-line prefer-destructuring
                output[pair[0]] = pair[1];
            });

            [numberArrayPairs, stringArrayPairs, booleanArrayPairs].forEach((arrayPairs) => {
                arrayPairs.forEach(([key, value]) => {
                    const [index, actualValue] = value;
                    if (!output[key]) output[key] = [];
                    output[key][index - 1] = actualValue; // to start arrays from index 0, while fortran lists start from 1
                });
            });

            return output;
        }

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
        // eslint-disable-next-line class-methods-use-this
        fortranExtractKeyValuePairs(data, regexPattern, type, isArray) {
            if (!typeParsers[type]) throw new Error("Invalid type");
            const parser = typeParsers[type];

            return Array.from(data.matchAll(regexPattern)).map((match) => {
                const key = match[1];
                const value = isArray
                    ? [parseInt(match[2], 10), parser(match[3])]
                    : parser(match[2]);

                return [key, value];
            });
        }
    };
