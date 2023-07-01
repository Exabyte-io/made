import s from "underscore.string";

// const fortranDoubleRegex = "([-+]?\\d*\\.?\\d*(?:[Eed][+-]?\\d+)?)";
// const fortranNamelistRegex = "&(%s)((?:\\s|\\S)*?)\\/"; // capturing group: 1 - namelist name, 2 - data
// const fortranCardsRegex = "^\\s*\\/(?![\\s\\S]*^\\/)([\\s\\S]*)"; // capturing group: 1 - cards section string
// const keyValueRegex = "^\\s*%s\\s*=\\s*%s\\s*\\n"; // 0 - string value name, 1 - str|bool|number
// const stringRegex = "([+\\w.\\-\\/]*)"; // capturing group: 1 - string value
// const fortranStringRegex = "'([\\w.\\-\\+\\/ ]*)'"; // capturing group: string value
// const fortranArrayRegex = "^\\s*%s\\(%s\\)\\s*=\\s*%s\\s*\\n"; // capturing groups: 0 - array name, 1 - array index, 2 - value
// const fortranBooleanRegex = "\\.(true|false)\\."; // capturing group: boolean value
const fortranDoubleRegex =
    "([-+]?" + // Optional leading sign
    "\\d*" + // Zero or more digits before the decimal point
    "\\.?" + // Optional decimal point
    "\\d*" + // Zero or more digits after the decimal point
    "(?:[Eed][+-]?\\d+)?" + // Optional exponent part
    ")";

const fortranNamelistRegex =
    "&" + // Start with an ampersand
    "%s" + // Namelist name
    "((?:\\s|\\S)*?)" + // Matches any sequence of space or non-space characters
    "\\/"; // Ends with a slash

const fortranCardsRegex =
    "^\\s*\\" + // A slash at the beginning of a line with any leading spaces
    "/" + // Negative lookahead
    "(?![\\s\\S]*^\\/)" + // Not followed by any characters then a line starting with a slash
    "([\\s\\S]*)"; // Capture all characters till end

const keyValueRegex =
    "^\\s*" + // Key name at the start of a line with any leading spaces
    "%s" + // Key name
    "\\s*=" + // Equal sign
    "\\s*" + // Value with any leading and trailing spaces
    "%s" + // Value
    "\\s*" +
    "\\n"; // Ends with a newline character

const stringRegex = "([+\\w.\\-\\/]*)"; // Matches alphanumeric, plus, period, hyphen, and slash characters

const fortranStringRegex =
    "'" + // Starting single quote
    "([\\w.\\-\\+\\/ ]*)" + // Matches alphanumeric, period, hyphen, plus, slash, and space characters
    "'"; // Ending single quote

const fortranArrayRegex =
    "^\\s*" + // Array name at the start of a line with any leading spaces
    "%s" + // Array name
    "\\(" + // Array index enclosed in parentheses
    "%s" + // Array index
    "\\)" +
    "\\s*=" + // Equal sign with any leading spaces
    "\\s*" + // Value with any leading and trailing spaces
    "%s" + // Value
    "\\s*" +
    "\\n"; // Ends with a newline character

const fortranBooleanRegex =
    "\\." + // Starting period
    "(true|false)" + // Matches either "true" or "false"
    "\\."; // Ending period

export const regex = {
    fortranDouble: fortranDoubleRegex,
    fortranNamelist: fortranNamelistRegex,
    fortranCards: fortranCardsRegex,
    keyValue: keyValueRegex,
    string: stringRegex,
    fortranString: fortranStringRegex,
    fortranArray: fortranArrayRegex,
    fortranBoolean: fortranBooleanRegex,

    stringKeyValue: new RegExp(s.sprintf(keyValueRegex, stringRegex, fortranStringRegex), "gm"),
    numberKeyValue: new RegExp(s.sprintf(keyValueRegex, stringRegex, fortranDoubleRegex), "gm"),
    booleanKeyValue: new RegExp(s.sprintf(keyValueRegex, stringRegex, fortranBooleanRegex), "gm"),
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
};
