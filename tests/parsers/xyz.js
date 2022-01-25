import { expect } from "chai";

import parsers from "../../src/parsers/parsers";
import { getNameFromContents } from "../../src/parsers/xyz";
import { GenericXYZ, Si } from "../enums";
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
    it("should return the material name based on the xyz file contents", () => {
        const fileContents = GenericXYZ;
        expect(getNameFromContents(fileContents)).to.be.equal("Methane");
    });
});
