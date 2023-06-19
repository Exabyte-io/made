import { formatString } from "../utils";

const fortranDoubleRegex = "([-+]?\\d*\\.?\\d*(?:[Eed][+-]?\\d+)?)";
const fortranNamelistRegex = "&({0})((?:\\s|\\S)*?)\\/"; // capturing group: 1 - namelist name, 2 - data
const fortranCardsRegex = "^\\s*\\/(?![\\s\\S]*^\\/)([\\s\\S]*)"; // capturing group: 1 - cards section string
const keyValueRegex = "^\\s*{0}\\s*=\\s*{1}\\s*\\n"; // 0 - string value name, 1 - str|bool|number
const stringRegex = "([+\\w.\\-\\/]*)"; // capturing group: 1 - string value
const fortranStringRegex = "'([\\w.\\-\\+\\/ ]*)'"; // capturing group: string value
const fortranArrayRegex = "^\\s*{0}\\({1}\\)\\s*=\\s*{2}\\s*\\n"; // capturing groups: 0 - array name, 1 - array index, 2 - value
const fortranBooleanRegex = "\\.(true|false)\\."; // capturing group: boolean value

export const regex = {
    fortranDouble: fortranDoubleRegex,
    fortranNamelist: fortranNamelistRegex,
    fortranCards: fortranCardsRegex,
    keyValue: keyValueRegex,
    string: stringRegex,
    fortranString: fortranStringRegex,
    fortranArray: fortranArrayRegex,
    fortranBoolean: fortranBooleanRegex,

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
};
