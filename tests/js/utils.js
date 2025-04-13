import { expect } from "chai";
import _ from "underscore";

export function assertDeepAlmostEqual(leftHandOperand, rightHandOperand, excludedKeys = []) {
    expect(_.omit(leftHandOperand, excludedKeys)).to.be.deep.almost.equal(
        _.omit(rightHandOperand, excludedKeys),
    );
}
