import { expect } from "chai";

import { Material } from "../../src/material";
import { ESPRESSOMaterialParser } from "../../src/parsers/espresso/parser";
import parsers from "../../src/parsers/parsers";
import { BNHex, BNHexIbravPWSCFInput, Si, SiPWSCFInput } from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Parsers:Espresso", () => {
    it("should return textual representation of a material according to QE pw.x input format", () => {
        const material = new Material(Si);
        expect(parsers.espresso.toEspressoFormat(material)).to.be.equal(SiPWSCFInput);
    });

    it("should return a material config from QE input file for BN", () => {
        const parser = new ESPRESSOMaterialParser();
        const materialConfig = parser.parse(BNHexIbravPWSCFInput);
        assertDeepAlmostEqual(materialConfig, BNHex, ["name"]);
    });
});
