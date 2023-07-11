import { expect } from "chai";

import { Material } from "../../src/material";
import { ESPRESSOMaterialParser } from "../../src/parsers/espresso/parser";
import parsers from "../../src/parsers/parsers";
import { BN, BNPWSCFInput, Si, SiPWSCFInput } from "../enums";

describe("Parsers:Espresso", () => {
    it("should return textual representation of a material according to QE pw.x input format", () => {
        const material = new Material(Si);
        expect(parsers.espresso.toEspressoFormat(material)).to.be.equal(SiPWSCFInput);
    });

    it("should return a material config from QE input file for BN", () => {
        const parser = new ESPRESSOMaterialParser();
        const materialConfig = parser.parse(BNPWSCFInput);
        console.log(materialConfig);
        expect(materialConfig).to.be.deep.equal(BN); // TODO: put actual material config from another commit
    });
});
