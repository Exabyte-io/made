import { expect } from "chai";

import parsers from "../../src/parsers/parsers";
import { validateAll } from "../../src/parsers/xyz";
import { CH4InvalidFormat, largeBasis, Si, smallBasis } from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Parsers:XYZ", () => {
    it("should extract basis from XYZ text", () => {
        const text = "Si 0 0 0 \n Si 0.25 0.25 0.25";
        assertDeepAlmostEqual(parsers.xyz.toBasisConfig(text), Si.basis, [
            "constraints",
            "cell",
            "units",
        ]);
    });
    it("should return error number for a structure having less than the maximum number of atoms", () => {
        const error = validateAll(smallBasis);
        expect(error).to.be.equal(0);
    });
    it("should return the error number for a structure having more than the maximum number of atoms", () => {
        const error = validateAll(largeBasis);
        expect(error).to.be.equal(2001);
    });
    it("should return the error number for a structure that does not have a valid format", () => {
        const error = validateAll(CH4InvalidFormat);
        expect(error).to.be.equal(1001);
    });
});
