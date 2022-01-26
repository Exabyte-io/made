import { expect } from "chai";

import parsers from "../../src/parsers/parsers";
import { validate } from "../../src/parsers/xyz";
import { largeBasis, Si, smallBasis } from "../enums";
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
        const error = validate(smallBasis);
        expect(error).to.be.equal(0);
    });
    it("should return the error number for a structure having more than the maximum number of atoms", () => {
        const error = validate(largeBasis);
        expect(error).to.be.equal(2001);
    });
});
