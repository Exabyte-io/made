import fs from "fs-extra";

export function readFile(filePath, coding = 'utf8') {
    return fs.readFileSync(filePath, coding);
}

export function readJSONFile(filePath) {
    return JSON.parse(readFile(filePath));
}
