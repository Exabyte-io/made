import s from "underscore.string";

const fortranDoubleRegex =
    "([-+]?" + // Optional leading sign
    "\\d*" + // Zero or more digits before the decimal point
    "\\.?" + // Optional decimal point
    "\\d*" + // Zero or more digits after the decimal point
    "(?:[EeDd][+-]?\\d+)?" + // Optional exponent part
    ")";

const fortranNamelistRegex =
    "&" + // Start with an ampersand
    "%s" + // Namelist name placeholder
    "((?:\\s|\\S)*?)" + // Matches any sequence of space or non-space characters
    "\\/"; // Ends with a slash

const fortranCardsRegex =
    "^\\s*\\/" + // Slash at the start of a line with any leading spaces
    "(?![\\s\\S]*^\\/)" + // Negative lookahead for a slash at the beginning of the next line
    "([\\s\\S]*)"; // Capture all characters till end

const keyValueRegex =
    "^\\s*" + // Key name at the start of a line with any leading spaces
    "%s" + // Key name placeholder
    "\\s*=\\s*" + // Equal sign with any leading and trailing spaces
    "%s" + // Value placeholder
    "\\s*\\n"; // Ends with a newline character

const fortranStringRegex =
    "'" + // Starting single quote
    "([\\w.\\-\\+\\/ ]*)" + // Matches alphanumeric, period, hyphen, plus, slash, and space characters
    "'"; // Ending single quote

const fortranArrayRegex =
    "^\\s*" + // Array name at the start of a line with any leading spaces
    "%s" + // Array name
    "\\(" + // Array index opening parentheses
    "%s" + // Array index
    "\\)" + // Array index closing parentheses
    "\\s*=\\s*" + // Equal sign with any leading and trailing spaces
    "%s" + // Value placeholder
    "\\s*\\n"; // Ends with a newline character

const fortranBooleanRegex =
    "\\." + // Starting period
    "(true|false)" + // Matches either "true" or "false" surrounded by periods
    "\\."; // Ending period

const stringRegex = "([+\\w.\\-\\/]*)"; // Matches alphanumeric, plus, period, hyphen, and slash characters
export const regex = {
    general: {
        double:
            "[-+]?" + // Optional leading sign
            "\\d*" + // Zero or more digits before the decimal point
            "\\.?" + // Optional decimal point
            "\\d*" + // Zero or more digits after the decimal point
            "(?:[Ee][+-]?\\d+)?", // Optional exponent part,
        string: stringRegex,
    },
    fortran: {
        stringKeyValue: new RegExp(s.sprintf(keyValueRegex, stringRegex, fortranStringRegex), "gm"),
        numberKeyValue: new RegExp(s.sprintf(keyValueRegex, stringRegex, fortranDoubleRegex), "gm"),
        booleanKeyValue: new RegExp(
            s.sprintf(keyValueRegex, stringRegex, fortranBooleanRegex),
            "gm",
        ),
        numberArrayKeyValue: new RegExp(
            s.sprintf(fortranArrayRegex, stringRegex, "(\\d+)", fortranDoubleRegex),
            "gm",
        ),
        stringArrayKeyValue: new RegExp(
            s.sprintf(fortranArrayRegex, stringRegex, "(\\d+)", fortranStringRegex),
            "gm",
        ),
        booleanArrayKeyValue: new RegExp(
            s.sprintf(fortranArrayRegex, stringRegex, "(\\d+)", fortranBooleanRegex),
            "gm",
        ),
        namelists: (namelistName) =>
            new RegExp(s.sprintf(fortranNamelistRegex, `(${namelistName.toUpperCase()})`)),
        cards: new RegExp(fortranCardsRegex, "m"),
    },
};
