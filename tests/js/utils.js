import { expect } from "chai";
import fs from "fs";
import _ from "underscore";

export function readFile(filePath, coding = "utf8") {
    return fs.readFileSync(filePath, coding);
}

export function readJSONFile(filePath) {
    return JSON.parse(readFile(filePath));
}

export function assertDeepAlmostEqual(leftHandOperand, rightHandOperand, excludedKeys = []) {
    expect(_.omit(leftHandOperand, excludedKeys)).to.be.deep.almost.equal(
        _.omit(rightHandOperand, excludedKeys),
    );
}
