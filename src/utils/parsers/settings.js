import s from "underscore.string";

const fortranDoubleRegex = "([-+]?\\d*\\.?\\d*(?:[Eed][+-]?\\d+)?)";
const fortranNamelistRegex = "&(%s)((?:\\s|\\S)*?)\\/"; // capturing group: 1 - namelist name, 2 - data
const fortranCardsRegex = "^\\s*\\/(?![\\s\\S]*^\\/)([\\s\\S]*)"; // capturing group: 1 - cards section string
const keyValueRegex = "^\\s*%s\\s*=\\s*%s\\s*\\n"; // 0 - string value name, 1 - str|bool|number
const stringRegex = "([+\\w.\\-\\/]*)"; // capturing group: 1 - string value
const fortranStringRegex = "'([\\w.\\-\\+\\/ ]*)'"; // capturing group: string value
const fortranArrayRegex = "^\\s*%s\\(%s\\)\\s*=\\s*%s\\s*\\n"; // capturing groups: 0 - array name, 1 - array index, 2 - value
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
        new RegExp(s.sprintf(fortranNamelistRegex, namelistName.toUpperCase())),
    cards: new RegExp(fortranCardsRegex, "m"),
};
