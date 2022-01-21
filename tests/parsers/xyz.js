import { expect } from "chai";

import parsers from "../../src/parsers/parsers";
import { validateNumberOfAtoms } from "../../src/parsers/xyz";
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
    it("should return the boolean of true for the file having less than the maximum number of atoms", () => {
        expect(validateNumberOfAtoms(smallBasis)).to.be.equal(true);
    });
    it("should return the boolean of false for the file having less than the maximum number of atoms", () => {
        expect(validateNumberOfAtoms(largeBasis)).to.be.equal(false);
    });
});
