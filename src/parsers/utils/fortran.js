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
                namelists[namelistName] = this.fortranExtractKeyValuePairs(namelistData);
            });
            return namelists;
        }

        /**
         * @summary Extracts an array of the key value pairs from a Fortran namelist.
         * @param {String} data - namelist data
         * @returns {Object}
         *
         * @example
         * for input data:
         *   ecutrho =   4.8000000000d+02
         *   ecutwfc =   6.0000000000d+01
         *   ibrav = 4
         *   celldm(1) = 4.7478008
         *   celldm(3) = 3.0676560682
         *   nat = 4
         *   nosym = .false.
         *   ntyp = 2
         *   occupations = 'fixed'
         *
         * should return object:
         *  {
         *    ecutrho: 480,
         *    ecutwfc: 60,
         *    ibrav: 4,
         *    celldm: [4.7478008, null, 3.0676560682],
         *    nat: 4,
         *    nosym: false,
         *    ntyp: 2,
         *    occupations: 'fixed'
         *    }
         */
        fortranExtractKeyValuePairs(data) {
            const pairTypes = [
                // TODO: add support for a number in form of 1.234D-56 -- current solution parses it as 1.234
                { regexPattern: regex.fortran.numberKeyValue, type: Number },
                { regexPattern: regex.fortran.stringKeyValue, type: String },
                { regexPattern: regex.fortran.booleanKeyValue, type: Boolean },
                // TODO: Change regex and capturing to accommodate for: Fortran lists assigned multiple values inline: list = 1,2,3 -- current implementation doesn't capture that
                { regexPattern: regex.fortran.numberArrayKeyValue, type: Number, isArray: true },
                { regexPattern: regex.fortran.stringArrayKeyValue, type: String, isArray: true },
                { regexPattern: regex.fortran.booleanArrayKeyValue, type: Boolean, isArray: true },
            ];

            return pairTypes.reduce((output, { regexPattern, type, isArray }) => {
                this.fortranExtractKeyValuePair(data, regexPattern, type, isArray).forEach(
                    ([key, value]) => {
                        if (isArray) {
                            output[key] = output[key] || [];
                            // eslint-disable-next-line prefer-destructuring
                            output[key][value[0] - 1] = value[1]; // to start arrays from index 0, while Fortran lists start from 1
                        } else {
                            output[key] = value;
                        }
                    },
                );
                return output;
            }, {});
        }

        /**
         * Extracts key-value pairs from a string data using provided regex pattern and type.
         * If isArray is set to true, treats the key-value pair as an array.
         *
         * @param {String} data - The string data to extract key-value pairs from.
         * @param {RegExp} regexPattern - The regex pattern to use for extracting.
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
        fortranExtractKeyValuePair(data, regexPattern, type, isArray = false) {
            const parser = typeParsers[type];
            return Array.from(data.matchAll(regexPattern)).map(([, key, index, value]) =>
                isArray ? [key, [parseInt(index, 10), parser(value)]] : [key, parser(index)],
            );
        }
    };
